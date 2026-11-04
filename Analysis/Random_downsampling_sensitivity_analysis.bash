$ python scripts/setup_sub.py DNA 80 --seed 1 --data_dir data/DNA_degapped --ref path_to_ref/HIV.HXB2.sequence.fasta 
---
A 199 -> 80
B 8743 -> 80
C 903 -> 80
D 89 -> 80
Rare 81 -> 80
CRF01.AE 1107 -> 80
CRF02.AG 92 -> 80
project written to proj_DNA_80_seed1
---

$ for s in $(seq 1 100); do   python scripts/setup_sub.py DNA 80 --seed $s --ref path_to_ref/HIV.HXB2.sequence.fasta; done
