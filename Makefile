.PHONY: pilot-large

# 40 events on 2 Mb through label.smk. Local only; run inside the container.
# Script archived under archive/scripts/ (B3); kept for optional dress rehearsal.
pilot-large:
	python archive/scripts/pilot_large.py
