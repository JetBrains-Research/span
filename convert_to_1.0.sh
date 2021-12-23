#!/bin/bash
set -eoup pipefail

if [ $# -lt 1 ];
then
  echo "Usage: ./convert_to_1.0.sh <old_span_model.span> [<new_span_model.span>]"
  exit 1
fi;

OLD_SPAN_MODEL=$1

if [ ! -f "$OLD_SPAN_MODEL" ];
then
  echo "File $OLD_SPAN_MODEL doesn't exist"
  exit 1
fi;

SPAN_TMPDIR=$(mktemp -d --tmpdir "span_converter_XXXXXX")

tar -C "$SPAN_TMPDIR" -xf "$OLD_SPAN_MODEL"

SPAN_TMP_INFO="$SPAN_TMPDIR/information.json"

if [ ! -f "$SPAN_TMP_INFO" ];
then
  echo "information.json missing"
  exit 1
fi;

# pre 0.13 version to 0.13
sed -E 's/      "path":/      "treatment":/g;s/  "version": [0-9]+/  "version": 3/;/  "bin_size": /a\  "fragment": "auto",\n  "fit.information.fqn": "org.jetbrains.bio.experiments.fit.Span1AnalyzeFitInformation",' -i "$SPAN_TMP_INFO"

# 0.13 to 1.0
sed -E 's/org.jetbrains.bio.experiments.fit/org.jetbrains.bio.span.fit/' -i "$SPAN_TMP_INFO"
sef -E 's/org.jetbrains.bio.statistics.hmm/org.jetbrains.bio.span.statistics.hmm' -i "$SPAN_TMP_INFO"

if [ $# -ge 2 ];
then
  NEW_SPAN_MODEL="$2"
else
  NEW_SPAN_MODEL="$OLD_SPAN_MODEL"
fi;

tar -cf "$NEW_SPAN_MODEL" -C "$SPAN_TMPDIR" .

rm -r "$SPAN_TMPDIR"

echo "Converted $OLD_SPAN_MODEL to $NEW_SPAN_MODEL"