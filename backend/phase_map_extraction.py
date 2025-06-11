import re
import csv

def extract_transitions_with_temps(phase_map_str: str):
    if not isinstance(phase_map_str, str):
        return []

    transitions = []
    segments = phase_map_str.split(";")

    valid_phases = {
        "CR", "N", "I", "FN", "NF", "SMA", "SMX", "FNG", "NX", "GLASS", "DECOMP", "SMZA"
    }

    for segment in segments:
        parts = re.split(r'\s*-\s*', segment.strip())
        parsed = []

        for part in parts:
            clean = part.strip().upper().rstrip("?")
            if clean in valid_phases:
                parsed.append(("phase", clean))
            elif re.match(r"[~>]?[-+]?[0-9]+(?:\.[0-9]+)?", clean):
                temp_match = re.search(r"[-+]?[0-9]+(?:\.[0-9]+)?", clean)
                if temp_match:
                    parsed.append(("temp", temp_match.group()))

        i = 0
        while i + 2 < len(parsed):
            if parsed[i][0] == "phase" and parsed[i+1][0] == "temp" and parsed[i+2][0] == "phase":
                phase1 = parsed[i][1]
                temp = parsed[i+1][1]
                phase2 = parsed[i+2][1]
                if phase1 != phase2:
                    transitions.append(f"{phase1} - {phase2} @ {temp}")
                i += 2
            else:
                i += 1

    return transitions


def parse_csv_to_transition_dict(csv_path):
    result_dict = {}
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            compound_id = row.get("id")
            phase_map = row.get("phase map", "")
            if compound_id:
                transitions = extract_transitions_with_temps(phase_map)
                result_dict[compound_id] = transitions
    return result_dict


if __name__ == "__main__":
    csv_path = "phase_map_edit_v2.csv"
    output_csv = "phase_map_output.csv"

    transitions_by_id = parse_csv_to_transition_dict(csv_path)

    # Show in terminal
    for compound_id, transitions in transitions_by_id.items():
        print(f"{compound_id}: {transitions}")

    # Write to output CSV
    with open(output_csv, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "parsed_transitions"])
        for cid, transitions in transitions_by_id.items():
            joined = "; ".join(transitions)
            writer.writerow([cid, joined])

    print(f"\n✅ Output written to {output_csv}")
