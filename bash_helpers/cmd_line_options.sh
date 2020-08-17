#!/bin/bash

for option; do
    echo $option
    case $option in
	-help | --help | -h)
	    want_help=yes
	    ;;

	-working-dir=* | --working-dir=*)
	    WORKINGDIR=`expr "x$option" : "x-*working-dir=\(.*\)"`
	    ;;

	# -wipe-existing-data=* | --wipe-existing-data=* )
	#     WIPEEXISTING=`expr "x$option" : "x-*wipe-existing-data=\(.*\)"`
	#     ;;

	# unrecognized option}
	-*)
	    { echo "error: unrecognized option: $option
	    	  Try \`$0 --help' for more information." >&2
	      { (exit 1); exit 1; }; }
	    ;;
    esac
done


if test "$want_help" = yes; then
  cat <<EOF
\`./main_tpls.sh' build tpls

Usage: $0 [OPTION]...

Options:
-h, --help			   display help and exit

--working-dir=			   the working directory where we build/run

EOF
  exit 0
fi

# --wipe-existing-data=[yes/no]	   if yes, all the following subfolders:
# 				     --target-dir/data_*
# 				     --target-dir/build_*
# 				   will be fully wiped.
# 				   default = no
