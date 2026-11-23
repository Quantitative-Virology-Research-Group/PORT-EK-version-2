#!/usr/bin/env python3
"""Run the ORIGINAL PORTEKrun.py with the in-memory corrections of portek_runtime_patch.py. Same arguments as PORTEKrun.py.
usage: python portek_patched_cli.py find_enriched PROJECT -k 15 -m 2.0   (set PORTEK_DIR to the repository, default ./PORT-EK-version-2)"""
import os, sys, runpy
repo = os.environ.get("PORTEK_DIR") or (os.getcwd() if os.path.exists("PORTEKrun.py") else "PORT-EK-version-2")
repo = os.path.abspath(repo)
if not os.path.exists(os.path.join(repo, "PORTEKrun.py")): sys.exit(f"PORTEKrun.py not found in {repo}; set PORTEK_DIR to the PORT-EK-version-2 folder"); sys.path.insert(0, repo); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import portek, portek_runtime_patch; portek_runtime_patch.apply_findk_patch()
sys.argv[0] = os.path.join(repo, "PORTEKrun.py"); runpy.run_path(sys.argv[0], run_name="__main__")
