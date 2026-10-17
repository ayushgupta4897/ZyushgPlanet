#!/usr/bin/env bash
# Fetch new Modal results, re-triage, and vet any survivor not yet vetted; loop until the batch is finished.
cd "$(dirname "$0")/.." && . .venv/bin/activate
RUN=${1:-poolA_main}
while true; do
  python -m pipeline.modal_client fetch --run-id "$RUN" 2>&1 | tail -1
  python -m pipeline.triage --run-id "$RUN" --config pipeline/prereg_poolA.json > /dev/null 2>&1
  python -m pipeline.vet_batch --run-id "$RUN" --skip-done 2>&1 | grep -E "STAGE2" | sed 's/.*| STAGE2/STAGE2/' | sort | uniq -c | sort -rn | head -5
  echo "$(date +%H:%M) vetted: $(ls runs/$RUN/vet/*/vet.json 2>/dev/null | wc -l) survivors; stage2 passes: $(grep -l '"stage2_pass": true' runs/$RUN/vet/*/vet.json 2>/dev/null | wc -l)"
  st=$(python -m pipeline.modal_client status --run-id "$RUN" 2>/dev/null | tr -d '\n ')
  echo "$st" | grep -q '"state":"finished"' && { echo "batch finished; final pass"; python -m pipeline.modal_client fetch --run-id "$RUN"; python -m pipeline.triage --run-id "$RUN" --config pipeline/prereg_poolA.json > /dev/null 2>&1; python -m pipeline.vet_batch --run-id "$RUN" --skip-done > /dev/null 2>&1; echo DONE; break; }
  sleep 600
done
