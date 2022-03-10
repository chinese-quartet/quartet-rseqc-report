#!/usr/bin/env bash

# A wrapper for exp2qcdt r package
#
# Author: Jingcheng Yang
# Email: yjcyxky@163.com
#
# License: MIT

# Exit on error. Append "|| true" if you expect an error.
# set -o errexit
# Exit on error inside any functions or subshells.
# set -o errtrace
# Do not allow use of undefined vars. Use ${VAR:-} to use an undefined VAR
set -o nounset
# Catch the error in case mysqldump fails (but gzip succeeds) in `mysqldump |gzip`
set -o pipefail
# Turn on traces, useful while debugging but commented out by default
# set -o xtrace

show_help() {
	cat <<EOF
usage: $(echo $0) [-e <FPKM_FILE>] [-c <COUNT_FILE>][-m <META_FILE>] [-o <RESULT_DIR>]
       -e FPKM_FILE Fpkm table file.
       -c COUNT_FILE Count table file
       -m META_FILE Metadata file.
       -o RESULT_DIR A directory for result files.
EOF
}

while getopts ":he:c:m:o:" arg; do
	case "$arg" in
	"e")
		FPKM_FILE="$OPTARG"
		;;
	"c")
		COUNT_FILE="$OPTARG"
		;;
	"m")
		META_FILE="$OPTARG"
		;;
	"o")
		RESULT_DIR="$OPTARG"
		;;
	"?")
		echo "Unkown option: $OPTARG"
		exit 1
		;;
	":")
		echo "No argument value for option $OPTARG"
		;;
	h)
		show_help
		exit 0
		;;
	*)
		echo "Unknown error while processing options"
		show_help
		exit 1
		;;
	esac
done

if [ -z "$FPKM_FILE" ]; then
	echo "-e argument is not specified."
	exit 1
elif [ ! -f "$FPKM_FILE" ]; then
	echo "$FPKM_FILE is not a valid file."
	exit 1
else
	FPKM_FILE=$(realpath $FPKM_FILE)
fi

if [ -z "$COUNT_FILE" ]; then
	echo "-c argument is not specified."
	exit 1
elif [ ! -f "$COUNT_FILE" ]; then
	echo "$COUNT_FILE is not a valid file."
	exit 1
else
	COUNT_FILE=$(realpath $COUNT_FILE)
fi

if [ -z "$META_FILE" ]; then
	echo "-m argument is not specified."
	exit 1
elif [ ! -f "$META_FILE" ]; then
	echo "$META_FILE is not a valid file."
	exit 1
else
	META_FILE=$(realpath $META_FILE)
fi

if [ -z "$RESULT_DIR" ]; then
	echo "-o argument is not specified."
	exit 1
elif [ ! -d "$RESULT_DIR" ]; then
	echo "$RESULT_DIR is not a valid directory."
	exit 1
else
	RESULT_DIR=$(realpath $RESULT_DIR)
fi

TEMP=$(mktemp)

echo "Run script: $TEMP"
cat <<EOF >"$TEMP"
#!/usr/bin/env Rscript

run <- function() {
	# Print traceback message
	on.exit(traceback())
	library(exp2qcdt)
	print("Running...")
	exp2qcdt("$FPKM_FILE", "$COUNT_FILE", "$META_FILE", "$RESULT_DIR")
}

run()
EOF

printf "\n---------------------\n"
cat "$TEMP"
echo "---------------------"

Rscript $TEMP
