#!/usr/bin/env python3
"""
Convert the MCParticles collection in an edm4hep2json-produced JSON file
into a Phoenix-native "Tracks" collection (straight line from vertex to
endpoint per particle), since Phoenix's EDM4hep-JSON importer does not
recognize edm4hep::MCParticleCollection as a visualizable type -- only
Tracks/Hits/Clusters/Jets/Vertices/etc. are natively supported.
(https://github.com/HSF/phoenix/blob/main/guides/developers/event_data_format.md)

Usage:
    ./mcparticles_to_phoenix_tracks.py input.edm4hep.json [-o output.json]

Produces a separate JSON file, structured the same "Event N" -> per-event
way as the input, with one top-level "Tracks" category per event holding
a single named collection ("MCParticleTracks"). Load it in Phoenix
alongside the original hits file -- both overlay in the same 3D view.
"""
import argparse
import json


def particle_to_track(p):
    v, e = p["vertex"], p["endpoint"]
    charge = p.get("charge", 0.0)
    mom = p.get("momentum", {"x": 0.0, "y": 0.0, "z": 0.0})
    mom_mag = (mom["x"] ** 2 + mom["y"] ** 2 + mom["z"] ** 2) ** 0.5
    # simple charge-based colouring: red = positive, blue = negative, grey = neutral
    color = "#ff4444" if charge > 0 else ("#4477ff" if charge < 0 else "#aaaaaa")
    return {
        "pos": [[v["x"], v["y"], v["z"]], [e["x"], e["y"], e["z"]]],
        "charge": charge,
        "mom": mom_mag,
        "color": color,
    }


def convert(infile, outfile):
    with open(infile) as f:
        data = json.load(f)

    out = {}
    for event_key, event in data.items():
        if not event_key.startswith("Event"):
            continue
        mcparticles = event.get("MCParticles", {}).get("collection", [])
        tracks = [particle_to_track(p) for p in mcparticles]
        out[event_key] = {"Tracks": {"MCParticleTracks": tracks}}

    with open(outfile, "w") as f:
        json.dump(out, f, indent=2)

    n_events = len(out)
    n_tracks = sum(len(v["Tracks"]["MCParticleTracks"]) for v in out.values())
    print(f"Wrote {outfile}: {n_events} event(s), {n_tracks} track(s) total")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("infile", help="edm4hep2json output file (must include MCParticles)")
    parser.add_argument("-o", "--out", default=None, help="output file (default: <infile>_tracks.json)")
    args = parser.parse_args()

    outfile = args.out or args.infile.rsplit(".", 1)[0] + "_tracks.json"
    convert(args.infile, outfile)
