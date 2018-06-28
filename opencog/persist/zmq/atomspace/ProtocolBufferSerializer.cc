/*
 * opencog/atomspace/ProtocolBufferSerializer.cc
 *
 * Copyright (C) 2008-2015 OpenCog Foundation
 * All Rights Reserved
 *
 * Written by Erwin Joosten, Hendy Irawan <ceefour666@gmail.com>
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

#include "ProtocolBufferSerializer.h"
#include "opencog/atoms/base/Handle.h"
#include "opencog/atoms/base/Atom.h"
#include "opencog/atoms/base/Link.h"
#include "opencog/atoms/base/Node.h"
#include "opencog/truthvalue/AttentionValue.h"
#include "opencog/truthvalue/DistributionalValue.h"
#include "opencog/atomspaceutils/TLB.h"

using namespace opencog;

TLB tlbuf;

ProtocolBufferSerializer::ProtocolBufferSerializer()
{
}

ProtocolBufferSerializer::~ProtocolBufferSerializer()
{
};


void ProtocolBufferSerializer::deserializeAtom(
        const ZMQAtomMessage& atomMessage, Atom& atom)
{
    //deserializeAttentionValueHolder(atomMessage.attentionvalueholder(), atom);

    // XXX FIXME This is truly a bad and illegal design -- screwing
    // around with the Atom internals is just plain the wrong thing to
    // do.  Please use the Atom ctors correctly!
    atom._atom_space = NULL;
    if (atomMessage.incoming_size() == 0)
    {
        atom._incoming_set = NULL;
    }
    else
    {
    	// TODO: is incoming set transferred over ZMQBackingStore? If so, how to "deserialize" them back?
    	atom._incoming_set = NULL;
//    	IncomingSet inSet(atomMessage.incoming_size());
//        atom._incoming_set = shared_ptr<LinkPtr>(inSet);
//        HandleEntry *previous = atom.incoming;
//        for(int i = 1; i < atomMessage.incoming_size(); i++)
//        {
//            HandleEntry *current = new HandleEntry(Handle(atomMessage.incoming(i)));
//        	LinkPtr linkPtr();
//        	inSet.push_back(linkPtr);
//        }
    }

    atom._flags = atomMessage.flags();
}

//void ProtocolBufferSerializer::serializeAtom(
//        Atom& atom, ZMQAtomMessage* atomMessage)
//{
//    //serializeAttentionValueHolder(atom, atomMessage->mutable_attentionvalueholder());
//
//    atomMessage->set_handle(atom.get_handle().value());
//
//    HandleEntry* next=atom.incoming;
//    while(next)
//    {
//        atomMessage->add_incoming(next->handle.value());
//        next = next->next;
//    }
//
//    atomMessage->set_type(atom.type);
//    atomMessage->set_flags(atom.flags);
//
//    serialize(*atom.truthValue, atomMessage->mutable_truthvalue());
//}

AtomPtr ProtocolBufferSerializer::deserialize(const ZMQAtomMessage& atomMessage)
{
    switch(atomMessage.atomtype())
    {
    case ZMQAtomTypeNode:
    {
    	NodePtr nodePtr = deserializeNode(atomMessage);
        return nodePtr;
    }
    case ZMQAtomTypeLink:
    {
        LinkPtr linkPtr = deserializeLink(atomMessage);
        return linkPtr;
    }
    case ZMQAtomTypeNotFound:
    	return AtomPtr();
    default:
        throw RuntimeException(TRACE_INFO, "Invalid ZMQ atomtype");
    }
}

//void ProtocolBufferSerializer::serialize(Atom &atom, ZMQAtomMessage* atomMessage)
//{
//    Link* link = dynamic_cast<Link *>(&atom);
//    if(link)
//        serializeLink(*link, atomMessage);
//    else
//    {
//        Node* node = dynamic_cast<Node *>(&atom);
//        if(node)
//            serializeNode(*node, atomMessage);
//        else
//            throw RuntimeException(TRACE_INFO, "Invalid atomtype");
//    }
//}

//void ProtocolBufferSerializer::deserializeAttentionValue(
//        const ZMQAttentionValueHolderMessage &attentionValueHolderMessage,
//        AttentionValue& av)
//{
//    av.m_STI=attentionValueHolderMessage.sti();
//    av.m_LTI=attentionValueHolderMessage.lti();
//    av.m_VLTI=attentionValueHolderMessage.vlti();
//}
//
//void ProtocolBufferSerializer::serializeAttentionValue(
//        AttentionValue& av, ZMQAttentionValueHolderMessage* attentionValueHolderMessage)
//{
//    attentionValueHolderMessage->set_sti(av.m_STI);
//    attentionValueHolderMessage->set_lti(av.m_LTI);
//    attentionValueHolderMessage->set_vlti(av.m_VLTI);
//}
//
//void ProtocolBufferSerializer::deserializeAttentionValueHolder(
//        const ZMQAttentionValueHolderMessage &attentionValueHolderMessage,
//        AttentionValueHolder& avh )
//{
//    deserializeAttentionValue(attentionValueHolderMessage, avh.attentionValue);
//}
//
//void ProtocolBufferSerializer::serializeAttentionValueHolder(
//        AttentionValueHolder& avh, ZMQAttentionValueHolderMessage* attentionValueHolderMessage)
//{
//    serializeAttentionValue(avh.attentionValue, attentionValueHolderMessage);
//}

NodePtr ProtocolBufferSerializer::deserializeNode(
        const ZMQAtomMessage& atomMessage)
{
	Handle nodePtr(createNode(atomMessage.type(), atomMessage.name()));
    if (atomMessage.has_truthvalue()) {
        DistributionalValuePtr tv;
    	tv = deserialize(atomMessage.truthvalue());
		nodePtr->setTruthValue(tv);
    }
	tlbuf.addAtom(nodePtr, atomMessage.handle());
    deserializeAtom(atomMessage, *nodePtr);

    return NodeCast(nodePtr);
}

LinkPtr ProtocolBufferSerializer::deserializeLink(
        const ZMQAtomMessage& atomMessage)
{
	HandleSeq oset(atomMessage.outgoing_size());
    for(int i = 0; i < atomMessage.outgoing_size(); i++)
    {
		// oset[i] = Handle(atomMessage.outgoing(i));
    }

	Handle linkPtr(createLink(oset, atomMessage.type()));
    if (atomMessage.has_truthvalue()) {
        DistributionalValuePtr tv;
    	tv = deserialize(atomMessage.truthvalue());
	    linkPtr->setTruthValue(tv);
    }
	tlbuf.addAtom(linkPtr, atomMessage.handle());
    deserializeAtom(atomMessage, *linkPtr);

    return LinkCast(linkPtr);
}

//void ProtocolBufferSerializer::serializeLink(
//        Link& link, ZMQAtomMessage * atomMessage)
//{
//    serializeAtom(link, atomMessage);
//
//    atomMessage->set_atomtype(ZMQAtomTypeLink);
//
//    for (Handle h : link.outgoing)
//        atomMessage->add_outgoing(h.value());
//
//    if(link.trail)
//        serializeTrail(*(link.trail), atomMessage->mutable_trail());
//}
//
//void ProtocolBufferSerializer::serializeNode(
//        Node& node, ZMQAtomMessage* atomMessage)
//{
//    serializeAtom(node, atomMessage);
//
//    atomMessage->set_atomtype(ZMQAtomTypeNode);
//
//    atomMessage->set_name(node.name);
//}

void ProtocolBufferSerializer::serialize(DistributionalValue &tv, ZMQTruthValueMessage* truthValueMessage)
{
    throw RuntimeException(TRACE_INFO, "Not Implemented.");
}

DistributionalValuePtr ProtocolBufferSerializer::deserialize(
        const ZMQTruthValueMessage& truthValueMessage)
{
    if (truthValueMessage.singletruthvalue_size() == 1)
        return deserialize(truthValueMessage.singletruthvalue(0));
    else
    {
    	throw std::runtime_error("CompositeTruthValue no longer supported");
    }
}

DistributionalValuePtr ProtocolBufferSerializer::deserialize(
        const ZMQSingleTruthValueMessage& singleTruthValueMessage)
{
    throw RuntimeException(TRACE_INFO, "Not Implemented.");
    switch (singleTruthValueMessage.truthvaluetype())
    {
         throw RuntimeException(TRACE_INFO, "Invalid ZMQ truthvaluetype: '%d'.",
                 singleTruthValueMessage.truthvaluetype());
    }
}
