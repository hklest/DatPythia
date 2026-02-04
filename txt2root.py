#!/usr/bin/env python
"""
Debug version of the script to parse Pythia final-state particles.

It prints out detailed debug information for each line:
  - line index
  - line content (trimmed)
  - token count
  - whether it's considered a 'header' line or a 'particle' line
  - whether it triggers the 'event finished' block

We store final-state particles (status==1). The first one with |PDG| == 11 is saved
as the electron; all others are the hadronic final state. Then we fill the TTree
when we see "Event finished".
"""

import sys
import ROOT
from array import array

def parse_pythia_final_state_debug(input_filename, output_filename):
    # Read the input file.
    with open(input_filename, 'r') as infile:
        lines = infile.readlines()

    # Create the output ROOT file and TTree.
    outfile = ROOT.TFile(output_filename, "RECREATE")
    tree = ROOT.TTree("events", "Final state particles (debug)")

    # --- Branches for the final state electron four–vector ---
    ele_px = array('f', [0.])
    ele_py = array('f', [0.])
    ele_pz = array('f', [0.])
    ele_E  = array('f', [0.])

    # --- Branches for the hadronic final state ---
    max_hadrons = 100
    n_hadrons = array('i', [0])
    had_pdgs = array('i', max_hadrons * [0])
    had_px   = array('f', max_hadrons * [0.])
    had_py   = array('f', max_hadrons * [0.])
    had_pz   = array('f', max_hadrons * [0.])
    had_E    = array('f', max_hadrons * [0.])

    # Create TTree branches
    tree.Branch("ele_px", ele_px, "ele_px/F")
    tree.Branch("ele_py", ele_py, "ele_py/F")
    tree.Branch("ele_pz", ele_pz, "ele_pz/F")
    tree.Branch("ele_E",  ele_E,  "ele_E/F")

    tree.Branch("n_hadrons", n_hadrons, "n_hadrons/I")
    tree.Branch("had_pdgs",  had_pdgs,  "had_pdgs[n_hadrons]/I")
    tree.Branch("had_px",    had_px,    "had_px[n_hadrons]/F")
    tree.Branch("had_py",    had_py,    "had_py[n_hadrons]/F")
    tree.Branch("had_pz",    had_pz,    "had_pz[n_hadrons]/F")
    tree.Branch("had_E",     had_E,     "had_E[n_hadrons]/F")

    # Keep track of how many events we fill
    total_events_filled = 0

    # Current event's final-state particles
    current_particles = []

    #print(f"Reading file: {input_filename}")
    #print(f"Number of lines = {len(lines)}")

    # We'll iterate with enumerate for better debugging
    for i, line in enumerate(lines):
        line_stripped = line.strip()
        tokens = line_stripped.split()
        n_tokens = len(tokens)

        # Debug output for each line
        #print(f"\nLine {i}: \"{line_stripped}\"")
        #print(f"  -> #tokens = {n_tokens}")

        if not line_stripped:
            #print("  -> Empty line, skipping.")
            continue

        # Check if line is a banner or separator
        if ("PYTHIA EVENT FILE" in line_stripped) or ("====" in line_stripped):
           # print("  -> It's a banner/separator line.")
            if "Event finished" in line_stripped:
                #print("  -> 'Event finished' found! Processing current event...")

                # Process the current event
                ele_px[0] = 0.0
                ele_py[0] = 0.0
                ele_pz[0] = 0.0
                ele_E[0]  = 0.0
                n_hadrons[0] = 0
                electron_found = False

                print(f"     #particles in current event: {len(current_particles)}")
                for p in current_particles:
                    pdg = p["pdg"]
                    if (abs(pdg) == 11) and (not electron_found):
                       # print(f"     Electron found with PDG={pdg}, px={p['px']}, py={p['py']}, pz={p['pz']}, E={p['E']}")
                        ele_px[0] = p["px"]
                        ele_py[0] = p["py"]
                        ele_pz[0] = p["pz"]
                        ele_E[0]  = p["E"]
                        electron_found = True
                    else:
                        if n_hadrons[0] < max_hadrons:
                            idx = n_hadrons[0]
                            had_pdgs[idx] = p["pdg"]
                            had_px[idx]   = p["px"]
                            had_py[idx]   = p["py"]
                            had_pz[idx]   = p["pz"]
                            had_E[idx]    = p["E"]
                            n_hadrons[0] += 1

                tree.Fill()
                total_events_filled += 1

                # Reset for next event
                current_particles = []
            else:
                #print("  -> Separator line but not 'Event finished'")
                continue

        # Optionally skip header lines if they have too many tokens
        # (Here we just guess header lines >= 25 tokens, but you can adjust.)
        header_token_threshold = 25
        if n_tokens >= header_token_threshold:
            #  print(f"  -> {n_tokens} tokens >= {header_token_threshold}; treating as header line. Skipping.")
            continue

        # Now try to parse as a particle line. We expect at least 10 tokens.
        if n_tokens < 10:
            # print("  -> <10 tokens, skipping.")
            continue

        # Attempt to parse the line as: token[1]: status, token[2]: PDG, tokens[6..9]: px, py, pz, E
        try:
            status = int(tokens[1])
            pdg    = int(tokens[2])
            px = float(tokens[6])
            py = float(tokens[7])
            pz = float(tokens[8])
            E  = float(tokens[9])
        except ValueError:
            #print("  -> Could not parse as integers/floats. Skipping.")
            continue

        # We only care about final-state (status==1)
        if status != 1:
            # print(f"  -> status = {status}, not final-state. Skipping.")
            continue

       # print(f"  -> final-state particle PDG={pdg}, px={px}, py={py}, pz={pz}, E={E}")
        current_particles.append({"pdg": pdg, "px": px, "py": py, "pz": pz, "E": E})

    # After the loop, check if there is a trailing event with no "Event finished"
    if current_particles:
       # print("\nReached end of file but still have un-finished event. Filling it now.")
        ele_px[0] = 0.0
        ele_py[0] = 0.0
        ele_pz[0] = 0.0
        ele_E[0]  = 0.0
        n_hadrons[0] = 0
        electron_found = False

        #print(f"  #particles in final event: {len(current_particles)}")
        for p in current_particles:
            #pdg = p["pdg"]
            if (abs(pdg) == 11) and (not electron_found):
                #print(f"   Electron found with PDG={pdg}, px={p['px']}, py={p['py']}, pz={p['pz']}, E={p['E']}")
                ele_px[0] = p["px"]
                ele_py[0] = p["py"]
                ele_pz[0] = p["pz"]
                ele_E[0]  = p["E"]
                electron_found = True
            else:
                if n_hadrons[0] < max_hadrons:
                    idx = n_hadrons[0]
                    had_pdgs[idx] = p["pdg"]
                    had_px[idx]   = p["px"]
                    had_py[idx]   = p["py"]
                    had_pz[idx]   = p["pz"]
                    had_E[idx]    = p["E"]
                    n_hadrons[0] += 1

        tree.Fill()
        total_events_filled += 1
        if(total_events_filled%10 == 0):
            print(f"Total Events Processed: {total_events_filled}")

    outfile.Write()
    outfile.Close()
    print(f"\nConversion complete. Output written to {output_filename}")
    print(f"Total events filled into TTree: {total_events_filled}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_pythia_to_root.py input.txt output.root")
        sys.exit(1)
    parse_pythia_final_state_debug(sys.argv[1], sys.argv[2])
