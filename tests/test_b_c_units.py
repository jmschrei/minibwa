"""Area B unit tests through C drivers: the ksw2 kernels against a brute-force
affine-gap DP, chaining against a brute-force chaining DP and hand-worked
anchors, hit processing and MAPQ with hand-computed values, the presets and
adaptive parameters, and the base-alignment bookkeeping in align.c and cs.c.
Each case's expected values and their derivation are in tests/c/drv_b_*.c."""

import pytest

import minibwa_helpers as H

KSW_CASES = [
	("global", "101", "2000", "0"),           # ksw_extz2_sse, affine
	("global", "102", "2000", "1"),           # ksw_extd2_sse, dual affine
	("global", "103", "800", "0", "1"),       # generic scoring matrix (transitions)
	("global", "104", "800", "1", "1"),
	("global", "120", "1000", "1", "0", "12", "3", "23", "1"),  # other dual-affine costs (long-gap threshold)
	("global", "121", "800", "1", "1", "6", "3", "30", "1"),
	("global", "122", "800", "0", "0", "6", "1", "6", "1"),
	("extension", "105", "2000", "0"),
	("extension", "106", "2000", "1"),
	("global_approx", "107", "1000", "0"),
	("global_approx", "108", "1000", "1"),
	("band", "109", "1500", "0", "2"),
	("band", "110", "1500", "1", "2"),
	("zdrop", "111", "0"),
	("zdrop", "112", "1"),
	("local", "113", "2000"),
	("gen_mat",),
	("push_cigar",),
	("extz2_limits",),
]


@pytest.mark.parametrize("case", KSW_CASES, ids=["-".join(c) for c in KSW_CASES])
def test_ksw2(case):
	"""ksw_extz2_sse/ksw_extd2_sse (global, extension, approximate max, banded, Z-drop),
	ksw_ll_i16/ksw_ll_u8_core and the ksw2.h helpers against independent references."""
	H.assert_c_case("drv_b_ksw", case[0], *case[1:])


CHAIN_CASES = [
	("lchain_hand",),
	("lchain_random", "201", "1500"),
	("lchain_random", "202", "1500"),
	("cal_mblen",),
	("gen_hit",),
	("split_hit",),
	("squeeze_a",),
	("sync_hits",),
]


@pytest.mark.parametrize("case", CHAIN_CASES, ids=["-".join(c) for c in CHAIN_CASES])
def test_chain(case):
	"""mb_lchain_dp on hand-worked anchors and against a brute-force chaining DP; the
	chain-to-hit conversions mb_cal_mblen, mb_gen_hit, mb_split_hit, mb_squeeze_a, mb_sync_hits."""
	H.assert_c_case("drv_b_chain", case[0], *case[1:])


HIT_CASES = ["set_parent", "set_sam_pri", "select_sub", "hit_sort", "filter_hits", "set_mapq", "high_cov",
	"presets", "opt_adap", "is_sr_mode", "helpers", "radix_sort"]


@pytest.mark.parametrize("case", HIT_CASES)
def test_hits(case):
	"""Primary/secondary assignment, secondary selection, sorting, filtering, MAPQ with
	hand-computed values, high-occurrence coverage, presets, adaptive parameters, hashes,
	mb_log2 and the radix sorts."""
	H.assert_c_case("drv_b_hits", case)


EXTRA_CASES = ["update_extra", "append_cigar", "write_tags", "update_dp_max"]


@pytest.mark.parametrize("case", EXTRA_CASES)
def test_extra(case):
	"""mb_update_extra (CIGAR clean-up, mlen/blen, segment score, --eqx, cs/ds/MD),
	mb_append_cigar, mb_write_cs_ds/mb_write_MD and mb_update_dp_max on hand-worked alignments."""
	H.assert_c_case("drv_b_extra", case)


def test_driver_case_lists():
	"""Every case a driver defines is run by this file (a new case cannot be forgotten)."""
	assert set(H.c_cases("drv_b_ksw")) == set(c[0] for c in KSW_CASES)
	assert set(H.c_cases("drv_b_chain")) == set(c[0] for c in CHAIN_CASES)
	assert set(H.c_cases("drv_b_hits")) == set(HIT_CASES)
	assert set(H.c_cases("drv_b_extra")) == set(EXTRA_CASES)
