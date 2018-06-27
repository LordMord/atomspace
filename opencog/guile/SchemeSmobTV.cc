/*
 * SchemeSmobTV.c
 *
 * Scheme small objects (SMOBS) for truth values.
 *
 * Copyright (c) 2008,2009 Linas Vepstas <linas@linas.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License v3 as
 * published by the Free Software Foundation and including the exceptions
 * at http://opencog.org/wiki/Licenses
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program; if not, write to:
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstddef>
#include <memory>
#include <libguile.h>

#include <opencog/truthvalue/DistributionalValue.h>
#include <opencog/guile/SchemeSmob.h>

using namespace opencog;

/* ============================================================== */

#ifdef USE_KEYWORD_LIST_NOT_USED
/**
 * Search for a truth value (demarked by #:tv) in a list of key-value
 * pairs.  Return the truth value if found, else return null.
 * Throw errors if the list is not strictly just key-value pairs
 *
 * XXX This code is not currently used, since it seems pointless
 * to have key-value pairs for this function. After all, an atom
 * can only have one truth value ever -- if we find a truth value, we
 * use it. We don't really need a key to tell us that its a truth value.
 * So punt, and get truth values implicitly. Meanwhile, this code is
 * stubbed out, for a rainy day, in case we need to resurrect key-value
 * pairs in the future.
 */
static DistributionalValuePtr get_tv_from_kvp(SCM kvp, const char * subrname, int pos)
{
	if (!scm_is_pair(kvp)) return NULL;

	do
	{
		SCM skey = SCM_CAR(kvp);

		// Verify that the first item is a keyword.
		if (!scm_is_keyword(skey))
			scm_wrong_type_arg_msg(subrname, pos, skey, "keyword");

		skey = scm_keyword_to_symbol(skey);
		skey = scm_symbol_to_string(skey);
		char * key = scm_to_utf8_string(skey);

		kvp = SCM_CDR(kvp);
		pos ++;
		if (!scm_is_pair(kvp))
		{
			scm_wrong_type_arg_msg(subrname, pos, kvp, "value following keyword");
		}

		if (0 == strcmp(key, "tv"))
		{
			SCM sval = SCM_CAR(kvp);
			scm_t_bits misctype = SCM_SMOB_FLAGS(sval);
			if (misctype != COG_SIMPLE_TV)
				scm_wrong_type_arg_msg(subrname, pos, sval, "opencog truth value");
			return scm_to_tv(sval);
		}
		free(key);

		kvp = SCM_CDR(kvp);
		pos ++;
	}
	while (scm_is_pair(kvp));

	return NULL;
}
#endif /* USE_KEYWORD_LIST_NOT_USED */

/**
 * Search for a truth value in a list of values.
 * Return the truth value if found, else return null.
 * Throw errors if the list is not strictly just key-value pairs
 */
DistributionalValuePtr SchemeSmob::get_tv_from_list(SCM slist)
{
	while (scm_is_pair(slist))
	{
		SCM sval = SCM_CAR(slist);
		if (SCM_SMOB_PREDICATE(SchemeSmob::cog_misc_tag, sval))
		{
			scm_t_bits misctype = SCM_SMOB_FLAGS(sval);
			switch (misctype)
			{
				case COG_PROTOM: {
					ProtoAtomPtr pa(scm_to_protom(sval));
					DistributionalValuePtr tv(TruthValueCast(pa));
					if (tv) return tv;
				}
				default:
					break;
			}
		}

		slist = SCM_CDR(slist);
	}

	return nullptr;
}

/* ============================================================== */

std::string SchemeSmob::tv_to_string(const DistributionalValuePtr& tv)
{
	return tv->to_string();
}

/* ============================================================== */

SCM SchemeSmob::tv_to_scm (const DistributionalValuePtr& tv)
{
	return protom_to_scm(ProtoAtomCast(tv));
}

/**
 * Create a new simple truth value, with indicated mean and confidence.
 *
SCM SchemeSmob::ss_new_stv (SCM smean, SCM sconfidence)
{
	double mean = scm_to_double(smean);
	double confidence = scm_to_double(sconfidence);

	DistributionalValuePtr tv = SimpleTruthValue::createTV(mean, confidence);
	return tv_to_scm(tv);
}
*/

/* ============================================================== */
/**
 * Return true if the scm is a truth value
 */
SCM SchemeSmob::ss_tv_p (SCM s)
{
	ProtoAtomPtr pa(scm_to_protom(s));
	if (nullptr == pa) return SCM_BOOL_F;

	if (classserver().isA(pa->get_type(), TRUTH_VALUE))
		return SCM_BOOL_T;

	scm_remember_upto_here_1(s);
	return SCM_BOOL_F;
}

/**
 * Return true if the scm is a truth value
 */
inline SCM SchemeSmob::tv_p (SCM s, Type wanted)
{
	ProtoAtomPtr pa(scm_to_protom(s));
	if (nullptr == pa) return SCM_BOOL_F;

	if (wanted == pa->get_type()) return SCM_BOOL_T;
	scm_remember_upto_here_1(s);
	return SCM_BOOL_F;
}

SCM SchemeSmob::ss_stv_p (SCM s)
{
	return tv_p(s, SIMPLE_TRUTH_VALUE);
}

SCM SchemeSmob::ss_ctv_p (SCM s)
{
	return tv_p(s, COUNT_TRUTH_VALUE);
}

SCM SchemeSmob::ss_itv_p (SCM s)
{
	return tv_p(s, INDEFINITE_TRUTH_VALUE);
}

SCM SchemeSmob::ss_ptv_p (SCM s)
{
	return tv_p(s, PROBABILISTIC_TRUTH_VALUE);
}

SCM SchemeSmob::ss_ftv_p (SCM s)
{
	return tv_p(s, FUZZY_TRUTH_VALUE);
}

SCM SchemeSmob::ss_etv_p (SCM s)
{
	return tv_p(s, EVIDENCE_COUNT_TRUTH_VALUE);
}

/* ============================================================== */

DistributionalValuePtr SchemeSmob::scm_to_tv(SCM tv)
{
	ProtoAtomPtr pa(scm_to_protom(tv));
	DistributionalValuePtr dv(TruthValueCast(pa));
	return dv;
}

DistributionalValuePtr SchemeSmob::verify_tv(SCM tv, const char *subrname, int pos)
{
	ProtoAtomPtr pa(scm_to_protom(tv));
	DistributionalValuePtr dv(TruthValueCast(pa));

	if (nullptr == dv)
		scm_wrong_type_arg_msg(subrname, pos, tv, "opencog truth value");

	return dv;
}

/**
 * Return association list holding contents of a truth value
 */
SCM SchemeSmob::ss_tv_get_value (SCM s)
{
	DistributionalValuePtr tv = verify_tv(s, "cog-tv->alist");
	Type tvt = tv->get_type();

	if (DISTRIBUTIONAL_VALUE == tvt)
	{
		SCM mode = scm_from_vector_double(tv->get_mode());
		SCM mean = scm_from_vector_double(tv->get_mean());
		SCM conf = scm_from_double(tv->get_confidence());
		SCM count = scm_from_double(tv->total_count());
		SCM smode = scm_from_utf8_symbol("mode");
		SCM smean = scm_from_utf8_symbol("mean");
		SCM sconf = scm_from_utf8_symbol("confidence");
		SCM scount = scm_from_utf8_symbol("count");

		SCM rc = SCM_EOL;
		rc = scm_acons(smode, mode, rc);
		rc = scm_acons(smean, mean, rc);
		rc = scm_acons(sconf, conf, rc);
		rc = scm_acons(scount, count, rc);
		scm_remember_upto_here_1(s);
		return rc;
	}
    if (CONDITIONAL_DISTRIBUTIONAL_VALUE == tvt)
    {
        throw RuntimeException(TRACE_INFO,"Not Implemented");
    }

	scm_remember_upto_here_1(s);
	return SCM_EOL;
}

/**
 * Return the truth value mean
 */
SCM SchemeSmob::ss_tv_get_mean(SCM s)
{
	DistributionalValuePtr tv = verify_tv(s, "cog-tv-mean");
	return scm_from_vector_double(tv->get_mean());
}

/**
 * Return the truth value mode
 */
SCM SchemeSmob::ss_tv_get_mode(SCM s)
{
	DistributionalValuePtr tv = verify_tv(s, "cog-tv-mode");
	return scm_from_vector_double(tv->get_mode());
}
/**
 * Return the truth value confidence
 */
SCM SchemeSmob::ss_tv_get_confidence(SCM s)
{
	DistributionalValuePtr tv = verify_tv(s, "cog-tv-confidence");
	return scm_from_double(tv->get_confidence());
}

/**
 * Return the truth value count
 */
SCM SchemeSmob::ss_tv_get_count(SCM s)
{
	DistributionalValuePtr tv = verify_tv(s, "cog-tv-count");
	return scm_from_double(tv->total_count());
}


SCM SchemeSmob::ss_tv_merge (SCM sta, SCM stb)
{
	DistributionalValuePtr tva = verify_tv(sta, "cog-tv-merge", 1);
	DistributionalValuePtr tvb = verify_tv(stb, "cog-tv-merge", 2);

	return tv_to_scm(tva->merge(tvb));
}

/* ===================== END OF FILE ============================ */
