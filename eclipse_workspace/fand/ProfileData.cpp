/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * ProfileData.cpp
 *
 *  Created on: Jul 22, 2010
 *      Author: gareth
 */

#include "ProfileData.h"
#include "util.h"

ProfileType
kitID(std::string kitname)
{
    to_upper(kitname);

    if(kitname == "IDENTIFILER")
    {
        return Identifiler;
    }
    else if (kitname == "POWERPLEX16")
    {
        return PowerPlex16;
    }
    else if (kitname == "SGM+")
    {
        return SGM_Plus;
    }
    else
    {
        return prof_type_unknown;
    }
}

std::string
kitName(ProfileType prof_type)
{
    switch(prof_type)
    {
        case Identifiler:
            return "Identifiler";

        case PowerPlex16:
            return "Powerplex16";

        case SGM_Plus:
            return "SGM+";

        case mixed:
            return "Mixed";

        case test:
            return "Test";

        case prof_type_unknown:
        default:
            return "";
    }
}

EvidenceType
evidenceID(std::string evidence_type)
{
    to_upper(evidence_type);

    if (evidence_type == "CRIME")
    {
        return crime;
    }
    else if (evidence_type == "REFERENCE" || evidence_type == "REF")
    {
        return reference;
    }
    else
    {
        return ev_type_unknown;
    }
}

std::string
evidenceName(EvidenceType ev_type)
{
    switch(ev_type)
    {
        case crime:
            return "CRIME";

        case reference:
            return "REF";

        case ev_type_unknown:
        default:
            return "";
    }
}

std::ostream &
operator<<(std::ostream &os, ProfileData const &p)
{
    os << "Kit=" << p.m_kit_type << " Evidence=" << p.m_evidence_type << " ID=" << p.m_id
       << " Dataset=" << p.m_dataset << " Sample ID=" << p.m_sample_id << " Profile ID=" << p.m_profile_id
       << " Num Contribs=" << p.m_num_contributors << " Error rate=" << p.m_error_rate << std::endl;

    return os;
}
