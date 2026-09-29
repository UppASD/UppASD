#!/bin/sh
# Run the maintainer-side pyswatter comparison on harness output.
#
# First run run_traj_checks.py with a real UppASD binary from a scratch
# directory.  It creates coord.oamtest.out and restart.oamtest.out under
# traj_checks/.  This script deliberately invokes pyswatter; it does not
# reimplement or approximate its observable.
set -eu

case_root=${1:-traj_checks}

for case_name in 'ell+1' 'shift_b' 'w_area'; do
    case_dir=$case_root/$case_name
    test -d "$case_dir"
    test -f "$case_dir/coord.oamtest.out"
    test -f "$case_dir/restart.oamtest.out"
    (
        cd "$case_dir"
        case "$case_name" in
            'ell+1'|'w_area')
                # UppASD's default origin for these centred cases is (20,20,0).
                shift_x=20.0
                shift_y=20.0
                ;;
            shift_b)
                # shift_b explicitly sets oam_origin to (0,0,0).  Its remaining
                # raw-frame difference is the documented C8 basis tilt; --m0 is
                # a ground-state component transform, not an origin correction.
                shift_x=0.0
                shift_y=0.0
                ;;
        esac
        pyswatter-animate spin-oam-balance coord.oamtest.out restart.oamtest.out \
            --lz-integration site --shift "$shift_x" "$shift_y" 0.0 --output ref.csv
    )
done
