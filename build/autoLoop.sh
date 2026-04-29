#!/bin/bash

for E in 2 5 10 15 20
do
cat > run_tmp.mac << EOF
/control/verbose 1
/run/verbose 1
/event/verbose 0
/tracking/verbose 0

/run/initialize
/analysis/setFileName ecal_${E}GeV
/gun/particle e-
/gun/energy ${E} GeV

/run/beamOn 100000
EOF

./ecalSim run_tmp.mac
done
