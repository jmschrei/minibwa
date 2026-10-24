"""Area D: helpers under the writers and the pipeline, through the C driver tests/c/drv_dkommon.c.

kom_sprintf_lite() writes every byte of SAM and PAF output and is compared
with snprintf(); kom_parse_num() reads numeric options; mb_escape() expands
-R and -H; l2b_meth_convert(), l2b_meth_rev() and l2b_intv2cid_meth()
implement the bisulfite index layout; kt_for() and kt_pipeline() run the
mapping and must process every item once and keep output order.
"""

import pytest

import minibwa_helpers as H


@pytest.mark.parametrize("case", ["sprintf_lite", "parse_num", "escape", "strdup", "meth_convert", "kt_for", "kt_pipeline"])
def test_kommon_case(case):
	"""kom_sprintf_lite/km_sprintf_lite match snprintf for %d %ld %u %s %c and literal text (returned length,
	NULL buffer, appends, growth); kom_parse_num reads k/m/g suffixes and decimals rounding to nearest and
	stops after the suffix; mb_escape expands \\t \\n \\r \\\\ in place; kom_strdup copies; C->T and G->A
	conversion of 0-3 codes and their strand reversal; kt_for processes each index exactly once with valid
	thread ids for 1-8 threads; kt_pipeline runs the steps of each item in order and outputs items in input
	order for 1-4 threads."""
	H.assert_c_case("drv_dkommon", case)


def test_intv2cid_meth(synth, index_of):
	"""l2b_intv2cid_meth maps an interval of the --meth BWT text (C->T forward, G->A forward, reverse
	complement of G->A, reverse complement of C->T) back to the contig, the forward-strand position, the
	strand and the conversion, at contig starts, middles and ends; intervals across a contig end or past the
	text return a negative value."""
	meth = index_of(synth.genome, meth=True)
	args = [meth, len(synth.genome.names)] + [synth.genome.length(n) for n in synth.genome.names]
	H.assert_c_case("drv_dkommon", "intv2cid_meth", *args)
