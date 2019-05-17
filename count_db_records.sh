WORKDIR=${1:-results}
for DB in ${WORKDIR}/*.db; do echo "__${DB}__"; litecli $DB -te "select (select count(*) from probes_seq_info) as probes_seq_info, (select count(*) from probes_filtered) as probes_filteredi, (select count(*) from probes_filtered where is_musicc=1) as filtered_musicc;"; done
