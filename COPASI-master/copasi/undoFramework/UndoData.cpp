/*
 * UndoData.cpp
 *
 *  Created on: 15 Aug 2014
 *      Author: dada
 */

#include "UndoData.h"

UndoData::UndoData() {
	// TODO Auto-generated constructor stub

}

UndoData::~UndoData() {
	// TODO Auto-generated destructor stub
}

std::string UndoData::getName() const
{
    return mName;
}

void UndoData::setName(std::string name)
{
    mName = name;
}

std::string UndoData::getType() const
{
    return mType;
}

void UndoData::setType(std::string type)
{
    mType = type;
}
