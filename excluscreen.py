#!/usr/bin/env python3
"""
DNA Primer Exclusivity Screener - Version 3
============================================
A Tkinter GUI application that analyzes DNA primer sequences using NCBI BLASTN,
identifies binding sites, and finds potential primer pair hits on the same sequences.

Transferred over to new project with version control: excluscreen_v2.

Author: Claude AI Assistant (Opus 4.5)
License: To be determined

Requirements:
    - Python 3.x
    - requests library (pip install requests)
    - tkinter (usually included with Python)

Usage:
    python excluscreen_3.py
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, filedialog, messagebox
import requests
import xml.etree.ElementTree as ET
import time
import csv
import re
import threading
import queue
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set
from urllib.parse import urlencode


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class Primer:
    """Represents a DNA primer sequence with its metadata."""
    name: str
    sequence: str
    index: int  # Position in input list (used for pairing)


@dataclass
class Hit:
    """
    Represents a primer hit found via BLAST.
    Simplified from v2's BindingSite - focuses on essential info only.
    """
    primer_name: str
    subject_id: str
    subject_accession: str
    subject_title: str
    strand: str  # 'plus' or 'minus'
    subject_start: int
    subject_end: int
    query_start: int
    query_end: int
    query_length: int
    alignment_length: int
    identity_count: int  # Number of identical bases
    identity_percent: float  # Percent identity
    coverage_percent: float  # Query coverage percentage
    mismatches: int
    gaps: int
    evalue: float
    bit_score: float
    taxid: Optional[int] = None


@dataclass
class HitPair:
    """
    Represents a pair of primer hits on the same subject sequence.
    This is the inclusive replacement for v2's Amplicon class.
    """
    subject_id: str
    subject_title: str
    primer1_name: str
    primer2_name: str
    distance: int  # Distance between the two hits
    primer1_strand: str
    primer2_strand: str
    primer1_start: int
    primer1_end: int
    primer2_start: int
    primer2_end: int
    primer1_identity: float
    primer2_identity: float
    primer1_coverage: float
    primer2_coverage: float
    primer1_mismatches: int
    primer2_mismatches: int
    orientation: str  # 'convergent', 'divergent', 'same_strand'
    taxid: Optional[int] = None


def get_normalized_pair_key(pair: HitPair) -> tuple:
    """
    Create a normalized key for hit pair deduplication.
    Ensures that pairs with the same primers (regardless of order) are treated as duplicates.
    """
    if pair.primer1_name <= pair.primer2_name:
        primer1, primer2 = pair.primer1_name, pair.primer2_name
        p1_id, p2_id = pair.primer1_identity, pair.primer2_identity
        p1_cov, p2_cov = pair.primer1_coverage, pair.primer2_coverage
        p1_mm, p2_mm = pair.primer1_mismatches, pair.primer2_mismatches
    else:
        primer1, primer2 = pair.primer2_name, pair.primer1_name
        p1_id, p2_id = pair.primer2_identity, pair.primer1_identity
        p1_cov, p2_cov = pair.primer2_coverage, pair.primer1_coverage
        p1_mm, p2_mm = pair.primer2_mismatches, pair.primer1_mismatches

    return (
        primer1, primer2,
        pair.distance,
        pair.orientation,
        round(p1_id, 1), round(p2_id, 1),
        round(p1_cov, 1), round(p2_cov, 1),
        p1_mm, p2_mm,
        pair.taxid
    )


# =============================================================================
# INPUT PARSING FUNCTIONS
# =============================================================================

def parse_primer_input(text: str) -> List[Primer]:
    """
    Parse primer sequences from user input.
    Supports both FASTA format and simple list format.
    """
    primers = []
    lines = text.strip().split('\n')

    # Check if FASTA format (starts with >)
    if any(line.strip().startswith('>') for line in lines):
        current_name = None
        current_seq = []

        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name and current_seq:
                    seq = ''.join(current_seq).upper()
                    if is_valid_dna(seq):
                        primers.append(Primer(
                            name=current_name,
                            sequence=seq,
                            index=len(primers)
                        ))
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.replace(' ', ''))

        if current_name and current_seq:
            seq = ''.join(current_seq).upper()
            if is_valid_dna(seq):
                primers.append(Primer(
                    name=current_name,
                    sequence=seq,
                    index=len(primers)
                ))
    else:
        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue

            parts = re.split(r'\s+', line, maxsplit=1)

            if len(parts) == 2 and is_valid_dna(parts[1]):
                primers.append(Primer(
                    name=parts[0],
                    sequence=parts[1].upper(),
                    index=len(primers)
                ))
            elif len(parts) == 1 and is_valid_dna(parts[0]):
                primers.append(Primer(
                    name=f"Primer_{len(primers) + 1}",
                    sequence=parts[0].upper(),
                    index=len(primers)
                ))

    return primers


def is_valid_dna(sequence: str) -> bool:
    """Check if a string is a valid DNA sequence."""
    valid_chars = set('ATCGRYSWKMBDHVN')
    return len(sequence) > 0 and all(c.upper() in valid_chars for c in sequence)


def parse_taxid_list(text: str) -> List[int]:
    """Parse a comma or newline-separated list of taxonomy IDs."""
    taxids = []
    parts = re.split(r'[,;\s]+', text.strip())

    for part in parts:
        part = part.strip()
        if part:
            try:
                taxid = int(part)
                if taxid > 0:
                    taxids.append(taxid)
            except ValueError:
                continue

    return taxids


def validate_email(email: str) -> bool:
    """Basic email validation."""
    pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    return bool(re.match(pattern, email.strip()))


# =============================================================================
# NCBI BLAST API FUNCTIONS
# =============================================================================

BLAST_PUT_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
BLAST_GET_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"


def get_user_agent(email: str) -> str:
    """Generate a User-Agent string with the user's email for NCBI compliance."""
    return f"PrimerExcluScreen/3.0 (Contact: {email})"


def build_entrez_query(excluded_taxids: List[int], refseq_only: bool = False) -> Optional[str]:
    """Build an ENTREZ_QUERY string to exclude specific taxonomy IDs and/or limit to RefSeq."""
    query_parts = []

    if refseq_only:
        query_parts.append("srcdb_refseq[PROP]")

    if excluded_taxids:
        exclusions = [f"txid{taxid}[ORGN]" for taxid in excluded_taxids]
        if len(exclusions) == 1:
            query_parts.append(f"NOT {exclusions[0]}")
        else:
            query_parts.append(f"NOT ({' OR '.join(exclusions)})")

    if not query_parts:
        return None

    return " AND ".join(query_parts)


def submit_blast_query(sequence: str, database: str = "nt",
                       email: str = "user@example.com",
                       excluded_taxids: List[int] = None,
                       refseq_only: bool = False,
                       word_size: int = 7,
                       evalue: float = 1000,
                       callback=None) -> Optional[str]:
    """
    Submit a BLASTN query to NCBI and return the Request ID (RID).

    V3 Changes: Simplified parameters, focus on inclusive search.
    """
    params = {
        'CMD': 'Put',
        'PROGRAM': 'blastn',
        'TASK': 'blastn-short',
        'DATABASE': database,
        'QUERY': sequence,
        'FORMAT_TYPE': 'XML',
        'HITLIST_SIZE': 5000,  # Increased for more inclusive results
        'WORD_SIZE': word_size,
        'EXPECT': evalue,
        'NUCL_PENALTY': -1,
        'NUCL_REWARD': 1,
        'GAPCOSTS': '5 2',
        'FILTER': 'F',  # No low-complexity filtering for short primers
    }

    entrez_query = build_entrez_query(excluded_taxids or [], refseq_only)
    if entrez_query:
        params['ENTREZ_QUERY'] = entrez_query
        if callback:
            callback(f"Applying filter: {entrez_query}")

    headers = {
        'User-Agent': get_user_agent(email)
    }

    try:
        if callback:
            callback(f"Submitting BLAST query for sequence ({len(sequence)} bp)...")

        response = requests.post(BLAST_PUT_URL, data=params, headers=headers)
        response.raise_for_status()

        text = response.text

        rid_match = re.search(r'RID = (\w+)', text)
        if rid_match:
            rid = rid_match.group(1)
            if callback:
                callback(f"BLAST query submitted. RID: {rid}")
            return rid
        else:
            if callback:
                callback("Error: Could not extract RID from BLAST response")
            return None

    except requests.RequestException as e:
        if callback:
            callback(f"Error submitting BLAST query: {str(e)}")
        return None


def poll_blast_results(rid: str, poll_interval: int = 30,
                       max_attempts: int = 60,
                       email: str = "user@example.com",
                       callback=None) -> Optional[str]:
    """Poll NCBI BLAST server for results until ready or timeout."""
    params = {
        'CMD': 'Get',
        'RID': rid,
        'FORMAT_TYPE': 'XML',
    }

    headers = {
        'User-Agent': get_user_agent(email)
    }

    for attempt in range(max_attempts):
        try:
            status_params = {
                'CMD': 'Get',
                'RID': rid,
                'FORMAT_OBJECT': 'SearchInfo',
            }

            status_response = requests.get(BLAST_GET_URL, params=status_params,
                                           headers=headers)
            status_response.raise_for_status()
            status_text = status_response.text

            if 'Status=WAITING' in status_text:
                if callback:
                    callback(f"BLAST search in progress... (attempt {attempt + 1}/{max_attempts})")
                time.sleep(poll_interval)
                continue
            elif 'Status=FAILED' in status_text:
                if callback:
                    callback("BLAST search failed on server side")
                return None
            elif 'Status=UNKNOWN' in status_text:
                if callback:
                    callback("BLAST RID expired or unknown")
                return None
            elif 'Status=READY' in status_text:
                if 'ThereAreHits=yes' in status_text or 'ThereAreHits=no' in status_text:
                    if callback:
                        callback("BLAST search complete. Fetching results...")

                    results_response = requests.get(BLAST_GET_URL, params=params,
                                                    headers=headers)
                    results_response.raise_for_status()
                    return results_response.text
            else:
                time.sleep(poll_interval)
                continue

        except requests.RequestException as e:
            if callback:
                callback(f"Network error polling BLAST: {str(e)}")
            time.sleep(poll_interval)
            continue

    if callback:
        callback("BLAST polling timed out")
    return None


def parse_blast_xml(xml_text: str, primer: Primer,
                    min_coverage: float = 50.0,
                    min_identity: float = 70.0,
                    excluded_taxids: Set[int] = None) -> List[Hit]:
    """
    Parse BLAST XML results and extract hits.

    V3 Changes:
    - Simplified filtering based on coverage and identity only
    - Removed 3' end matching logic
    - More inclusive hit collection
    """
    hits = []
    excluded_taxids = excluded_taxids or set()

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        print(f"XML parsing error: {e}")
        return hits

    iterations = root.findall('.//Iteration')

    for iteration in iterations:
        blast_hits = iteration.findall('.//Hit')

        for blast_hit in blast_hits:
            hit_id = blast_hit.find('Hit_id')
            hit_def = blast_hit.find('Hit_def')
            hit_accession = blast_hit.find('Hit_accession')

            subject_id = hit_accession.text if hit_accession is not None else (
                hit_id.text if hit_id is not None else "Unknown")

            subject_accession = hit_accession.text if hit_accession is not None else "Unknown"

            subject_title = ""
            if hit_def is not None and hit_def.text:
                subject_title = hit_def.text

            # Try to extract taxid (best effort)
            taxid = extract_taxid_from_hit(hit_id, hit_def)

            if taxid and taxid in excluded_taxids:
                continue

            hsps = blast_hit.findall('.//Hsp')

            for hsp in hsps:
                try:
                    query_from = int(hsp.find('Hsp_query-from').text)
                    query_to = int(hsp.find('Hsp_query-to').text)
                    hit_from = int(hsp.find('Hsp_hit-from').text)
                    hit_to = int(hsp.find('Hsp_hit-to').text)

                    qseq = hsp.find('Hsp_qseq').text
                    hseq = hsp.find('Hsp_hseq').text

                    # Calculate statistics
                    alignment_length = len(qseq)

                    # Count identities
                    identity_count = sum(1 for q, h in zip(qseq, hseq)
                                         if q.upper() == h.upper() and q != '-' and h != '-')

                    # Calculate mismatches and gaps
                    mismatches = sum(1 for q, h in zip(qseq, hseq)
                                     if q != h and q != '-' and h != '-')
                    gaps = qseq.count('-') + hseq.count('-')

                    # Calculate percentages
                    identity_percent = (identity_count / alignment_length * 100) if alignment_length > 0 else 0

                    # Coverage based on query length
                    query_aligned_length = abs(query_to - query_from) + 1
                    coverage_percent = (query_aligned_length / len(primer.sequence) * 100)

                    # Get e-value and bit score
                    evalue_elem = hsp.find('Hsp_evalue')
                    evalue = float(evalue_elem.text) if evalue_elem is not None else 999

                    bitscore_elem = hsp.find('Hsp_bit-score')
                    bit_score = float(bitscore_elem.text) if bitscore_elem is not None else 0

                    # Determine strand
                    if hit_from <= hit_to:
                        strand = 'plus'
                    else:
                        strand = 'minus'

                    # Apply filtering thresholds
                    if coverage_percent >= min_coverage and identity_percent >= min_identity:
                        hits.append(Hit(
                            primer_name=primer.name,
                            subject_id=subject_id,
                            subject_accession=subject_accession,
                            subject_title=subject_title[:100] if len(subject_title) > 100 else subject_title,
                            strand=strand,
                            subject_start=min(hit_from, hit_to),
                            subject_end=max(hit_from, hit_to),
                            query_start=query_from,
                            query_end=query_to,
                            query_length=len(primer.sequence),
                            alignment_length=alignment_length,
                            identity_count=identity_count,
                            identity_percent=round(identity_percent, 1),
                            coverage_percent=round(coverage_percent, 1),
                            mismatches=mismatches,
                            gaps=gaps,
                            evalue=evalue,
                            bit_score=bit_score,
                            taxid=taxid
                        ))
                except (AttributeError, ValueError, TypeError) as e:
                    # Skip malformed HSPs
                    continue

    return hits


def extract_taxid_from_hit(hit_id, hit_def) -> Optional[int]:
    """Attempt to extract taxonomy ID from BLAST hit information."""
    # This is a placeholder - standard BLAST XML doesn't always include taxid directly
    return None


# =============================================================================
# HIT PAIR DETECTION FUNCTIONS
# =============================================================================

def determine_orientation(hit1: Hit, hit2: Hit) -> str:
    """
    Determine the relative orientation of two primer hits.

    Returns:
        'convergent': primers face each other (potential amplicon)
        'divergent': primers face away from each other
        'same_strand': both on same strand
    """
    if hit1.strand == hit2.strand:
        return 'same_strand'

    # Get the midpoints of each hit
    mid1 = (hit1.subject_start + hit1.subject_end) / 2
    mid2 = (hit2.subject_start + hit2.subject_end) / 2

    if hit1.strand == 'plus' and hit2.strand == 'minus':
        if mid1 < mid2:
            return 'convergent'  # Plus strand hit is upstream, facing downstream
        else:
            return 'divergent'
    else:  # hit1 minus, hit2 plus
        if mid2 < mid1:
            return 'convergent'
        else:
            return 'divergent'


def calculate_distance(hit1: Hit, hit2: Hit) -> int:
    """
    Calculate the distance between two hits.
    Returns the gap between the hits (0 if overlapping).
    """
    # Get the range of each hit
    range1 = (hit1.subject_start, hit1.subject_end)
    range2 = (hit2.subject_start, hit2.subject_end)

    # Check for overlap
    if range1[1] >= range2[0] and range2[1] >= range1[0]:
        return 0  # Overlapping

    # Calculate gap
    if range1[1] < range2[0]:
        return range2[0] - range1[1]
    else:
        return range1[0] - range2[1]


def find_hit_pairs(hits1: List[Hit], hits2: List[Hit],
                   max_distance: int,
                   same_primer: bool = False) -> List[HitPair]:
    """
    Find pairs of primer hits on the same subject within the specified distance.

    This is the V3 replacement for find_amplicons - much more inclusive.
    It finds ANY two hits on the same entry regardless of strand orientation.

    Args:
        hits1: Hits from first primer
        hits2: Hits from second primer
        max_distance: Maximum distance between hits
        same_primer: True if hits1 and hits2 are from the same primer

    Returns:
        List of HitPair objects
    """
    pairs = []

    # Group hits by subject ID
    hits1_by_subject: Dict[str, List[Hit]] = {}
    hits2_by_subject: Dict[str, List[Hit]] = {}

    for hit in hits1:
        if hit.subject_id not in hits1_by_subject:
            hits1_by_subject[hit.subject_id] = []
        hits1_by_subject[hit.subject_id].append(hit)

    for hit in hits2:
        if hit.subject_id not in hits2_by_subject:
            hits2_by_subject[hit.subject_id] = []
        hits2_by_subject[hit.subject_id].append(hit)

    # Find common subjects
    common_subjects = set(hits1_by_subject.keys()) & set(hits2_by_subject.keys())

    for subject in common_subjects:
        subject_hits1 = hits1_by_subject[subject]
        subject_hits2 = hits2_by_subject[subject]

        for i, h1 in enumerate(subject_hits1):
            # For same primer pairs, only compare with hits after this one to avoid duplicates
            start_j = i + 1 if same_primer else 0

            for j, h2 in enumerate(subject_hits2[start_j:] if same_primer else subject_hits2):
                # Skip if comparing identical hit (same position, same strand)
                if same_primer and h1.subject_start == h2.subject_start and h1.subject_end == h2.subject_end:
                    continue

                distance = calculate_distance(h1, h2)

                if distance <= max_distance:
                    orientation = determine_orientation(h1, h2)

                    pairs.append(HitPair(
                        subject_id=subject,
                        subject_title=h1.subject_title or h2.subject_title,
                        primer1_name=h1.primer_name,
                        primer2_name=h2.primer_name,
                        distance=distance,
                        primer1_strand=h1.strand,
                        primer2_strand=h2.strand,
                        primer1_start=h1.subject_start,
                        primer1_end=h1.subject_end,
                        primer2_start=h2.subject_start,
                        primer2_end=h2.subject_end,
                        primer1_identity=h1.identity_percent,
                        primer2_identity=h2.identity_percent,
                        primer1_coverage=h1.coverage_percent,
                        primer2_coverage=h2.coverage_percent,
                        primer1_mismatches=h1.mismatches,
                        primer2_mismatches=h2.mismatches,
                        orientation=orientation,
                        taxid=h1.taxid or h2.taxid
                    ))

    # Sort by distance
    pairs.sort(key=lambda x: x.distance)

    return pairs


# =============================================================================
# CSV EXPORT FUNCTION
# =============================================================================

def export_to_csv(pairs: List[HitPair], filepath: str) -> bool:
    """Export hit pair results to a CSV file."""
    try:
        # Deduplicate pairs for export
        unique_pairs = {}
        for pair in pairs:
            key = get_normalized_pair_key(pair)
            if key not in unique_pairs:
                unique_pairs[key] = {'example': pair.subject_id, 'title': pair.subject_title, 'count': 1}
            else:
                unique_pairs[key]['count'] += 1
                # Update title if current one is empty and we have a better one
                if not unique_pairs[key].get('title') and pair.subject_title:
                    unique_pairs[key]['title'] = pair.subject_title

        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)

            writer.writerow([
                'Example_Accession',
                'Accession_Name',
                'Accession_Count',
                'Primer_1',
                'Primer_2',
                'Distance',
                'Orientation',
                'P1_Identity_%',
                'P2_Identity_%',
                'P1_Coverage_%',
                'P2_Coverage_%',
                'P1_Mismatches',
                'P2_Mismatches'
            ])

            for key, info in unique_pairs.items():
                primer1, primer2, distance, orientation, p1_id, p2_id, p1_cov, p2_cov, p1_mm, p2_mm, taxid = key
                writer.writerow([
                    info['example'],
                    info.get('title', 'N/A')[:100] if info.get('title') else 'N/A',
                    info['count'],
                    primer1,
                    primer2,
                    distance,
                    orientation,
                    p1_id,
                    p2_id,
                    p1_cov,
                    p2_cov,
                    p1_mm,
                    p2_mm
                ])

        return True
    except IOError as e:
        print(f"Error writing CSV: {e}")
        return False


# =============================================================================
# MAIN GUI APPLICATION CLASS
# =============================================================================

class PrimerBlastApp:
    """
    Main application class for the DNA Primer Exclusivity Screener V3.
    """

    def __init__(self, root: tk.Tk):
        """Initialize the application."""
        self.root = root
        self.root.title("DNA Primer Exclusivity Screener v3")
        self.root.geometry("1000x850")
        self.root.minsize(900, 750)

        # Store results
        self.all_pairs: List[HitPair] = []
        self.all_hits: Dict[str, List[Hit]] = {}

        # Analysis state
        self.is_running = False

        # Message queue for thread-safe GUI updates
        self.gui_queue = queue.Queue()

        # Build the GUI
        self.create_widgets()

        # Start processing GUI update queue
        self._process_gui_queue()

    def create_widgets(self):
        """Create and arrange all GUI widgets."""

        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky="nsew")

        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(1, weight=1)
        main_frame.rowconfigure(4, weight=2)

        # === Primer Input Section ===
        input_header_frame = ttk.Frame(main_frame)
        input_header_frame.grid(row=0, column=0, sticky="ew", pady=(0, 5))
        input_header_frame.columnconfigure(0, weight=1)
        
        input_label = ttk.Label(input_header_frame,
                                text="Primer Sequences (FASTA format or one per line):",
                                font=('TkDefaultFont', 10, 'bold'))
        input_label.grid(row=0, column=0, sticky="w")
        
        self.import_button = ttk.Button(input_header_frame, text="Import FASTA File",
                                        command=self.import_fasta_file)
        self.import_button.grid(row=0, column=1, sticky="e", padx=(10, 0))

        self.primer_text = scrolledtext.ScrolledText(main_frame, height=8, width=80,
                                                     wrap=tk.WORD)
        self.primer_text.grid(row=1, column=0, sticky="nsew", pady=(0, 10))

        placeholder = """>Primer_1
ATCGATCGATCGATCGATCG
>Primer_2
TAGCTAGCTAGCTAGCTAGC
>Primer_3
GCGCGCGCATATATATAT
>Primer_4
ATATATATGCGCGCGCGC"""
        self.primer_text.insert("1.0", placeholder)
        self.primer_text.config(fg='gray')

        self.primer_text.bind('<FocusIn>', self.on_primer_focus_in)
        self.primer_text.bind('<FocusOut>', self.on_primer_focus_out)

        # === Parameters Section ===
        params_frame = ttk.LabelFrame(main_frame, text="Analysis Parameters (V3 - Inclusive Search)",
                                      padding="10")
        params_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))

        for i in range(6):
            params_frame.columnconfigure(i, weight=1 if i % 2 == 1 else 0)

        # Row 0: Hit detection parameters
        ttk.Label(params_frame, text="Min Coverage %:").grid(
            row=0, column=0, sticky="e", padx=(0, 5))
        self.min_coverage_var = tk.StringVar(value="50")
        self.min_coverage_entry = ttk.Entry(params_frame, textvariable=self.min_coverage_var,
                                            width=10)
        self.min_coverage_entry.grid(row=0, column=1, sticky="w", padx=(0, 20))

        ttk.Label(params_frame, text="Min Identity %:").grid(
            row=0, column=2, sticky="e", padx=(0, 5))
        self.min_identity_var = tk.StringVar(value="70")
        self.min_identity_entry = ttk.Entry(params_frame, textvariable=self.min_identity_var,
                                            width=10)
        self.min_identity_entry.grid(row=0, column=3, sticky="w", padx=(0, 20))

        ttk.Label(params_frame, text="Max Distance (bp):").grid(
            row=0, column=4, sticky="e", padx=(0, 5))
        self.max_distance_var = tk.StringVar(value="5000")
        self.max_distance_entry = ttk.Entry(params_frame, textvariable=self.max_distance_var,
                                            width=10)
        self.max_distance_entry.grid(row=0, column=5, sticky="w")

        # Row 1: BLAST database and email
        ttk.Label(params_frame, text="BLAST Database:").grid(
            row=1, column=0, sticky="e", padx=(0, 5), pady=(10, 0))
        self.db_var = tk.StringVar(value="nt")
        self.db_combo = ttk.Combobox(params_frame, textvariable=self.db_var,
                                     values=["nt", "nr", "refseq_genomic", "refseq_rna"],
                                     state="readonly", width=15)
        self.db_combo.grid(row=1, column=1, sticky="w", pady=(10, 0))

        ttk.Label(params_frame, text="Email (required):").grid(
            row=1, column=2, sticky="e", padx=(0, 5), pady=(10, 0))
        self.email_var = tk.StringVar(value="")
        self.email_entry = ttk.Entry(params_frame, textvariable=self.email_var,
                                     width=30)
        self.email_entry.grid(row=1, column=3, columnspan=3, sticky="w", pady=(10, 0))

        # Row 2: Taxonomy exclusion
        ttk.Label(params_frame, text="Exclude TaxIDs:").grid(
            row=2, column=0, sticky="ne", padx=(0, 5), pady=(10, 0))

        taxid_frame = ttk.Frame(params_frame)
        taxid_frame.grid(row=2, column=1, columnspan=5, sticky="w", pady=(10, 0))

        self.taxid_var = tk.StringVar(value="")
        self.taxid_entry = ttk.Entry(taxid_frame, textvariable=self.taxid_var,
                                     width=40)
        self.taxid_entry.pack(side="left")

        taxid_help = ttk.Label(taxid_frame,
                               text="  (comma-separated, e.g., 9606 for human, 10090 for mouse)",
                               font=('TkDefaultFont', 8), foreground='gray')
        taxid_help.pack(side="left", padx=(5, 0))

        # Row 3: RefSeq filter checkbox
        self.refseq_only_var = tk.BooleanVar(value=False)
        self.refseq_checkbox = ttk.Checkbutton(
            params_frame,
            text="Limit to RefSeq sequences only",
            variable=self.refseq_only_var
        )
        self.refseq_checkbox.grid(row=3, column=0, columnspan=3, sticky="w", pady=(10, 0))

        refseq_help = ttk.Label(params_frame,
                                text="(filters results to curated NCBI Reference Sequences)",
                                font=('TkDefaultFont', 8), foreground='gray')
        refseq_help.grid(row=3, column=3, columnspan=3, sticky="w", pady=(10, 0))

        # === Control Buttons Section ===
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=3, column=0, sticky="ew", pady=(0, 10))

        self.run_button = ttk.Button(button_frame, text="Run Analysis",
                                     command=self.run_analysis)
        self.run_button.pack(side="left", padx=(0, 10))

        self.export_button = ttk.Button(button_frame, text="Export to CSV",
                                        command=self.export_results,
                                        state="disabled")
        self.export_button.pack(side="left", padx=(0, 10))

        self.clear_button = ttk.Button(button_frame, text="Clear Results",
                                       command=self.clear_results)
        self.clear_button.pack(side="left")

        # Progress bar
        self.progress_var = tk.DoubleVar(value=0)
        self.progress_bar = ttk.Progressbar(button_frame, variable=self.progress_var,
                                            maximum=100, length=200)
        self.progress_bar.pack(side="right", padx=(10, 0))

        # === Results Section ===
        results_label = ttk.Label(main_frame, text="Results:",
                                  font=('TkDefaultFont', 10, 'bold'))
        results_label.grid(row=4, column=0, sticky="nw", pady=(0, 5))

        self.results_notebook = ttk.Notebook(main_frame)
        self.results_notebook.grid(row=4, column=0, sticky="nsew", pady=(25, 10))

        # Hit Pairs tab (replaces Amplicons)
        pairs_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(pairs_frame, text="Hit Pairs")
        pairs_frame.columnconfigure(0, weight=1)
        pairs_frame.rowconfigure(1, weight=1)

        # Filter frame for convergent pairs only
        filter_frame = ttk.Frame(pairs_frame)
        filter_frame.grid(row=0, column=0, sticky="ew", pady=(0, 5))
        
        self.filter_convergent_var = tk.BooleanVar(value=False)
        self.filter_checkbox = ttk.Checkbutton(
            filter_frame,
            text="Show only convergent pairs",
            variable=self.filter_convergent_var,
            command=self.update_displayed_pairs
        )
        self.filter_checkbox.pack(side="left")
        
        # Treeview frame with scrollbars
        pairs_tree_frame = ttk.Frame(pairs_frame)
        pairs_tree_frame.grid(row=1, column=0, sticky="nsew")
        pairs_tree_frame.columnconfigure(0, weight=1)
        pairs_tree_frame.rowconfigure(0, weight=1)
        
        columns = ('subject', 'subject_name', 'count', 'primer1', 'primer2', 'distance',
                   'orientation', 'p1_id', 'p2_id', 'p1_cov', 'p2_cov')
        self.pairs_tree = ttk.Treeview(pairs_tree_frame, columns=columns, show='headings',
                                       selectmode='browse')

        self.pairs_tree.heading('subject', text='Example Accession')
        self.pairs_tree.heading('subject_name', text='Accession Name')
        self.pairs_tree.heading('count', text='Count')
        self.pairs_tree.heading('primer1', text='Primer 1')
        self.pairs_tree.heading('primer2', text='Primer 2')
        self.pairs_tree.heading('distance', text='Distance')
        self.pairs_tree.heading('orientation', text='Orientation')
        self.pairs_tree.heading('p1_id', text='P1 Id%')
        self.pairs_tree.heading('p2_id', text='P2 Id%')
        self.pairs_tree.heading('p1_cov', text='P1 Cov%')
        self.pairs_tree.heading('p2_cov', text='P2 Cov%')

        self.pairs_tree.column('subject', width=120, minwidth=80)
        self.pairs_tree.column('subject_name', width=200, minwidth=100)
        self.pairs_tree.column('count', width=50, minwidth=40)
        self.pairs_tree.column('primer1', width=80, minwidth=60)
        self.pairs_tree.column('primer2', width=80, minwidth=60)
        self.pairs_tree.column('distance', width=70, minwidth=50)
        self.pairs_tree.column('orientation', width=80, minwidth=60)
        self.pairs_tree.column('p1_id', width=60, minwidth=45)
        self.pairs_tree.column('p2_id', width=60, minwidth=45)
        self.pairs_tree.column('p1_cov', width=60, minwidth=45)
        self.pairs_tree.column('p2_cov', width=60, minwidth=45)

        pairs_yscroll = ttk.Scrollbar(pairs_tree_frame, orient="vertical",
                                      command=self.pairs_tree.yview)
        pairs_xscroll = ttk.Scrollbar(pairs_tree_frame, orient="horizontal",
                                      command=self.pairs_tree.xview)
        self.pairs_tree.configure(yscrollcommand=pairs_yscroll.set,
                                  xscrollcommand=pairs_xscroll.set)
        
        self.pairs_tree.grid(row=0, column=0, sticky="nsew")
        pairs_yscroll.grid(row=0, column=1, sticky="ns")
        pairs_xscroll.grid(row=1, column=0, sticky="ew")

        # Status/Log tab
        log_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(log_frame, text="Status Log")
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        self.status_text = scrolledtext.ScrolledText(log_frame, height=15,
                                                     state='disabled')
        self.status_text.grid(row=0, column=0, sticky="nsew")

        # Individual Hits tab (replaces Binding Sites)
        hits_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(hits_frame, text="Individual Hits")
        hits_frame.columnconfigure(0, weight=1)
        hits_frame.rowconfigure(0, weight=1)

        self.hits_text = scrolledtext.ScrolledText(hits_frame, height=15,
                                                   state='disabled')
        self.hits_text.grid(row=0, column=0, sticky="nsew")

    def on_primer_focus_in(self, event):
        """Handle focus in on primer text - clear placeholder."""
        if self.primer_text.get("1.0", "end-1c").startswith(">Primer_1"):
            self.primer_text.delete("1.0", tk.END)
            self.primer_text.config(fg='black')

    def on_primer_focus_out(self, event):
        """Handle focus out on primer text - restore placeholder if empty."""
        if not self.primer_text.get("1.0", "end-1c").strip():
            placeholder = """>Primer_1
ATCGATCGATCGATCGATCG
>Primer_2
TAGCTAGCTAGCTAGCTAGC"""
            self.primer_text.insert("1.0", placeholder)
            self.primer_text.config(fg='gray')
    
    def import_fasta_file(self):
        """Import primer sequences from a FASTA file."""
        filepath = filedialog.askopenfilename(
            title="Import FASTA File",
            filetypes=[("FASTA files", "*.fa *.fasta *.fas *.txt"), ("All files", "*.*")]
        )
        
        if filepath:
            try:
                with open(filepath, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                if content.strip():
                    self.primer_text.config(state='normal', fg='black')
                    self.primer_text.delete("1.0", tk.END)
                    self.primer_text.insert("1.0", content)
                    self.primer_text.config(state='normal')
                    self.log_status(f"Imported FASTA file: {filepath}")
                else:
                    messagebox.showwarning("Empty File", "The selected file is empty.")
            except Exception as e:
                messagebox.showerror("Import Error", f"Failed to import file:\n{str(e)}")

    def log_status(self, message: str):
        """Add a message to the status log (thread-safe)."""
        if threading.current_thread() is threading.main_thread():
            self._do_log_status(message)
        else:
            self.gui_queue.put(('log_status', message))

    def _do_log_status(self, message: str):
        """Actually perform the log status update (must be called from main thread)."""
        self.status_text.config(state='normal')
        timestamp = time.strftime("%H:%M:%S")
        self.status_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.status_text.see(tk.END)
        self.status_text.config(state='disabled')

    def _process_gui_queue(self):
        """Process pending GUI updates from the queue."""
        try:
            while True:
                action, *args = self.gui_queue.get_nowait()

                if action == 'log_status':
                    self._do_log_status(args[0])
                elif action == 'set_progress':
                    self.progress_var.set(args[0])
                elif action == 'display_hits':
                    self._do_display_hits(args[0], args[1])
                elif action == 'display_pairs':
                    self._do_display_pairs(args[0])
                elif action == 'analysis_complete':
                    self._on_analysis_complete(args[0], args[1], args[2])
                elif action == 'analysis_error':
                    self._on_analysis_error(args[0])

        except queue.Empty:
            pass

        self.root.after(100, self._process_gui_queue)

    def _set_inputs_enabled(self, enabled: bool):
        """Enable or disable all input controls."""
        state = 'normal' if enabled else 'disabled'

        self.primer_text.config(state=state)
        self.min_coverage_entry.config(state=state)
        self.min_identity_entry.config(state=state)
        self.max_distance_entry.config(state=state)
        self.email_entry.config(state=state)
        self.taxid_entry.config(state=state)
        self.db_combo.config(state=state)
        self.refseq_checkbox.config(state=state)
        self.run_button.config(state=state)
        self.clear_button.config(state=state)

    def validate_parameters(self) -> Tuple[bool, float, float, int, str, List[int], bool]:
        """Validate user input parameters."""
        try:
            min_coverage = float(self.min_coverage_var.get())
            min_identity = float(self.min_identity_var.get())
            max_distance = int(self.max_distance_var.get())

            if min_coverage < 0 or min_coverage > 100:
                raise ValueError("Coverage must be between 0 and 100")
            if min_identity < 0 or min_identity > 100:
                raise ValueError("Identity must be between 0 and 100")
            if max_distance <= 0:
                raise ValueError("Max distance must be positive")

        except ValueError as e:
            messagebox.showerror("Invalid Parameters",
                                 f"Please enter valid numeric parameters.\n{str(e)}")
            return False, 0, 0, 0, "", [], False

        email = self.email_var.get().strip()
        if not email:
            messagebox.showerror("Email Required",
                                 "Please enter your email address.\n\n"
                                 "NCBI requires a valid email for BLAST queries.")
            return False, 0, 0, 0, "", [], False

        if not validate_email(email):
            messagebox.showerror("Invalid Email",
                                 "Please enter a valid email address.")
            return False, 0, 0, 0, "", [], False

        taxid_text = self.taxid_var.get().strip()
        excluded_taxids = parse_taxid_list(taxid_text) if taxid_text else []

        refseq_only = self.refseq_only_var.get()

        return True, min_coverage, min_identity, max_distance, email, excluded_taxids, refseq_only

    def run_analysis(self):
        """Main analysis workflow - triggered by Run Analysis button."""

        if self.is_running:
            messagebox.showinfo("Analysis Running",
                                "An analysis is already in progress.")
            return

        self.clear_results()

        valid, min_coverage, min_identity, max_distance, email, excluded_taxids, refseq_only = self.validate_parameters()
        if not valid:
            return

        primer_text = self.primer_text.get("1.0", tk.END)
        if primer_text.strip().startswith(">Primer_1"):
            messagebox.showwarning("No Input",
                                   "Please enter your primer sequences.")
            return

        primers = parse_primer_input(primer_text)

        if len(primers) < 1:
            messagebox.showwarning("No Primers",
                                   "Please enter at least 1 primer.")
            return

        num_pair_combinations = len(primers) * (len(primers) - 1) // 2
        num_self_combinations = len(primers)
        total_combinations = num_pair_combinations + num_self_combinations
        self.log_status(
            f"Parsed {len(primers)} primers ({num_pair_combinations} unique pairs + {num_self_combinations} self-tests = {total_combinations} combinations)")
        self.log_status(f"Using email: {email}")
        self.log_status(
            f"Parameters: Coverage >= {min_coverage}%, Identity >= {min_identity}%, Max Distance = {max_distance} bp")

        if excluded_taxids:
            self.log_status(f"Excluding taxonomy IDs: {', '.join(map(str, excluded_taxids))}")

        if refseq_only:
            self.log_status("RefSeq filter: ENABLED")

        database = self.db_var.get()

        self.is_running = True
        self._set_inputs_enabled(False)

        analysis_thread = threading.Thread(
            target=self._run_analysis_thread,
            args=(primers, database, email, excluded_taxids, refseq_only,
                  min_coverage, min_identity, max_distance),
            daemon=True
        )
        analysis_thread.start()

    def _run_analysis_thread(self, primers: List[Primer], database: str, email: str,
                             excluded_taxids: List[int], refseq_only: bool,
                             min_coverage: float, min_identity: float, max_distance: int):
        """Worker thread that performs the actual BLAST analysis."""
        try:
            hits_by_primer: Dict[str, List[Hit]] = {}
            excluded_taxids_set = set(excluded_taxids)

            total_primers = len(primers)

            for i, primer in enumerate(primers):
                self.log_status(f"\n--- Processing primer {i + 1}/{total_primers}: {primer.name} ---")
                self.log_status(
                    f"    Sequence: {primer.sequence[:30]}{'...' if len(primer.sequence) > 30 else ''} ({len(primer.sequence)} bp)")
                self.gui_queue.put(('set_progress', (i / total_primers) * 100))

                rid = submit_blast_query(
                    primer.sequence,
                    database,
                    email=email,
                    excluded_taxids=excluded_taxids,
                    refseq_only=refseq_only,
                    callback=self.log_status
                )

                if not rid:
                    self.log_status(f"Failed to submit BLAST for {primer.name}")
                    continue

                if i < total_primers - 1:
                    self.log_status("Waiting 2 seconds before next submission (NCBI courtesy)...")
                    time.sleep(2)

                xml_results = poll_blast_results(
                    rid,
                    poll_interval=30,
                    email=email,
                    callback=self.log_status
                )

                if not xml_results:
                    self.log_status(f"Failed to retrieve results for {primer.name}")
                    continue

                hits = parse_blast_xml(
                    xml_results,
                    primer,
                    min_coverage=min_coverage,
                    min_identity=min_identity,
                    excluded_taxids=excluded_taxids_set
                )
                hits_by_primer[primer.name] = hits

                self.log_status(
                    f"Found {len(hits)} hits for {primer.name} (coverage >= {min_coverage}%, identity >= {min_identity}%)")

                self.gui_queue.put(('display_hits', primer.name, hits))

            # Find hit pairs for ALL primer combinations
            self.log_status("\n=== Finding Hit Pairs (All Primer Combinations) ===")
            all_pairs = []

            num_combinations = len(primers) * (len(primers) - 1) // 2
            combination_count = 0

            for i in range(len(primers)):
                for j in range(i + 1, len(primers)):
                    primer1 = primers[i]
                    primer2 = primers[j]
                    combination_count += 1

                    self.log_status(
                        f"\nTesting combination {combination_count}/{num_combinations}: {primer1.name} vs {primer2.name}")

                    hits1 = hits_by_primer.get(primer1.name, [])
                    hits2 = hits_by_primer.get(primer2.name, [])

                    if not hits1:
                        self.log_status(f"  No hits for {primer1.name}")
                        continue
                    if not hits2:
                        self.log_status(f"  No hits for {primer2.name}")
                        continue

                    pairs = find_hit_pairs(hits1, hits2, max_distance)

                    if pairs:
                        self.log_status(f"  Found {len(pairs)} hit pairs within {max_distance} bp")
                        all_pairs.extend(pairs)
                    else:
                        self.log_status(f"  No pairs found within {max_distance} bp")

            # Test each primer against itself
            self.log_status("\n=== Testing Self-Pairing (Same Primer Hits) ===")
            for primer in primers:
                hits = hits_by_primer.get(primer.name, [])

                if len(hits) < 2:
                    self.log_status(f"  {primer.name}: Fewer than 2 hits, skipping self-test")
                    continue

                self.log_status(f"\nTesting self-pairing: {primer.name}")

                self_pairs = find_hit_pairs(hits, hits, max_distance, same_primer=True)

                if self_pairs:
                    self.log_status(f"  Found {len(self_pairs)} self-pairs")
                    all_pairs.extend(self_pairs)
                else:
                    self.log_status(f"  No self-pairs found")

            self.gui_queue.put(('analysis_complete', all_pairs, hits_by_primer, None))

        except Exception as e:
            import traceback
            self.gui_queue.put(('analysis_error', f"{str(e)}\n{traceback.format_exc()}"))

    def _on_analysis_complete(self, all_pairs: List[HitPair],
                              hits_by_primer: Dict[str, List[Hit]],
                              _unused):
        """Handle analysis completion."""
        self.all_pairs = all_pairs
        self.all_hits = hits_by_primer

        self._do_display_pairs(all_pairs)  # This will now use update_displayed_pairs internally

        # Calculate unique patterns
        unique_patterns = set()
        for pair in all_pairs:
            key = get_normalized_pair_key(pair)
            unique_patterns.add(key)

        self.progress_var.set(100)
        self._do_log_status(f"\n=== Analysis Complete ===")
        self._do_log_status(f"Total hit pairs found: {len(all_pairs)}")
        self._do_log_status(f"Unique patterns (after deduplication): {len(unique_patterns)}")

        self.is_running = False
        self._set_inputs_enabled(True)

        if all_pairs:
            self.export_button.config(state='normal')
            messagebox.showinfo("Analysis Complete",
                                f"Found {len(all_pairs)} hit pairs.\n"
                                f"({len(unique_patterns)} unique patterns after deduplication)\n"
                                f"See Results tab for details.")
        else:
            messagebox.showinfo("Analysis Complete",
                                "No hit pairs found within the specified parameters.")

    def _on_analysis_error(self, error_message: str):
        """Handle analysis error."""
        self._do_log_status(f"Error during analysis: {error_message}")

        self.is_running = False
        self._set_inputs_enabled(True)

        messagebox.showerror("Analysis Error", error_message)

    def display_hits(self, primer_name: str, hits: List[Hit]):
        """Display hits (thread-safe)."""
        if threading.current_thread() is threading.main_thread():
            self._do_display_hits(primer_name, hits)
        else:
            self.gui_queue.put(('display_hits', primer_name, hits))

    def _do_display_hits(self, primer_name: str, hits: List[Hit]):
        """Actually display hits (must be called from main thread)."""
        self.hits_text.config(state='normal')

        # Deduplicate hits by key parameters
        unique_hits = {}
        for hit in hits:
            key = (
                hit.strand,
                hit.subject_start,
                hit.subject_end,
                round(hit.identity_percent, 1),
                round(hit.coverage_percent, 1),
                hit.mismatches
            )
            if key not in unique_hits:
                unique_hits[key] = {'example': hit.subject_id, 'title': hit.subject_title, 'count': 1}
            else:
                unique_hits[key]['count'] += 1

        self.hits_text.insert(tk.END,
                              f"\n=== {primer_name} ({len(hits)} total hits, {len(unique_hits)} unique patterns) ===\n")

        display_count = 0
        for key, info in list(unique_hits.items())[:30]:
            strand, subj_start, subj_end, identity, coverage, mismatches = key
            example_accession = info['example']
            title = info['title']
            count = info['count']

            count_str = f" ({count} accessions)" if count > 1 else ""
            title_display = title[:40] + "..." if len(title) > 40 else title

            self.hits_text.insert(tk.END,
                                  f"  {example_accession}{count_str}\n"
                                  f"    {title_display}\n"
                                  f"    Strand: {strand}, Pos: {subj_start}-{subj_end}\n"
                                  f"    Identity: {identity}%, Coverage: {coverage}%, Mismatches: {mismatches}\n"
                                  )
            display_count += 1

        if len(unique_hits) > 30:
            self.hits_text.insert(tk.END, f"  ... and {len(unique_hits) - 30} more unique patterns\n")

        self.hits_text.see(tk.END)
        self.hits_text.config(state='disabled')

    def display_pairs(self, pairs: List[HitPair]):
        """Display hit pairs (thread-safe)."""
        if threading.current_thread() is threading.main_thread():
            self._do_display_pairs(pairs)
        else:
            self.gui_queue.put(('display_pairs', pairs))

    def _do_display_pairs(self, pairs: List[HitPair]):
        """Actually display hit pairs (must be called from main thread)."""
        # Store all pairs for filtering
        self.all_pairs = pairs
        self.update_displayed_pairs()
    
    def update_displayed_pairs(self):
        """Update the displayed pairs based on the filter setting."""
        for item in self.pairs_tree.get_children():
            self.pairs_tree.delete(item)

        # Apply filter if enabled
        pairs_to_display = self.all_pairs
        if self.filter_convergent_var.get():
            pairs_to_display = [p for p in self.all_pairs if p.orientation == 'convergent']

        unique_pairs = {}
        for pair in pairs_to_display:
            key = get_normalized_pair_key(pair)
            if key not in unique_pairs:
                unique_pairs[key] = {'example': pair.subject_id, 'title': pair.subject_title, 'count': 1}
            else:
                unique_pairs[key]['count'] += 1

        for key, info in unique_pairs.items():
            primer1, primer2, distance, orientation, p1_id, p2_id, p1_cov, p2_cov, p1_mm, p2_mm, taxid = key
            example_accession = info['example']
            subject_title = info.get('title', 'N/A')
            count = info['count']

            subject_display = example_accession[:20] + "..." if len(example_accession) > 20 else example_accession
            title_display = subject_title[:50] + "..." if len(subject_title) > 50 else (subject_title if subject_title else 'N/A')

            self.pairs_tree.insert('', tk.END, values=(
                subject_display,
                title_display,
                count,
                primer1,
                primer2,
                distance,
                orientation,
                p1_id,
                p2_id,
                p1_cov,
                p2_cov
            ))

    def export_results(self):
        """Export hit pair results to CSV file (exports what is currently displayed)."""
        # Get currently displayed pairs (respecting filter)
        pairs_to_export = self.all_pairs
        if self.filter_convergent_var.get():
            pairs_to_export = [p for p in self.all_pairs if p.orientation == 'convergent']
        
        if not pairs_to_export:
            messagebox.showwarning("No Results", "No hit pairs to export.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            title="Export Results"
        )

        if filepath:
            success = export_to_csv(pairs_to_export, filepath)
            if success:
                filter_status = " (filtered: convergent pairs only)" if self.filter_convergent_var.get() else ""
                self.log_status(f"Results exported to: {filepath}{filter_status}")
                messagebox.showinfo("Export Complete",
                                    f"Results exported to:\n{filepath}\n\n"
                                    f"Exported {len(pairs_to_export)} hit pairs.")
            else:
                messagebox.showerror("Export Failed",
                                     "Failed to export results. Check file permissions.")

    def clear_results(self):
        """Clear all results and reset the interface."""
        for item in self.pairs_tree.get_children():
            self.pairs_tree.delete(item)

        self.status_text.config(state='normal')
        self.status_text.delete("1.0", tk.END)
        self.status_text.config(state='disabled')

        self.hits_text.config(state='normal')
        self.hits_text.delete("1.0", tk.END)
        self.hits_text.config(state='disabled')

        self.progress_var.set(0)

        self.all_pairs = []
        self.all_hits = {}

        self.export_button.config(state='disabled')


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    """Main entry point for the application."""
    root = tk.Tk()

    try:
        style = ttk.Style()
        available_themes = style.theme_names()
        if 'clam' in available_themes:
            style.theme_use('clam')
        elif 'vista' in available_themes:
            style.theme_use('vista')
    except Exception:
        pass

    app = PrimerBlastApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()