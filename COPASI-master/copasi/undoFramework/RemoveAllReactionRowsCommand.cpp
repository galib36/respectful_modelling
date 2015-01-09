/*
 * RemoveAllReactionRowsCommand.cpp
 *
 *  Created on: 12 Aug 2014
 *      Author: dada
 */


#include <QtCore/QList>

#include "report/CCopasiRootContainer.h"
#include "model/CReaction.h"
#include "model/CReactionInterface.h"
#include "model/CModel.h"
#include "CQReactionDM.h"


#include "RemoveAllReactionRowsCommand.h"
#include "UndoReactionData.h"


RemoveAllReactionRowsCommand::RemoveAllReactionRowsCommand(CQReactionDM * pReaDM, const QModelIndex&) {
	mpReactionDM = pReaDM;

	assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
	CCopasiDataModel* pDataModel = (*CCopasiRootContainer::getDatamodelList())[0];
	assert(pDataModel != NULL);
	CModel * pModel = pDataModel->getModel();

	assert(pModel != NULL);

	for (int i = 0; i != pReaDM->rowCount()-1; ++i)
	{
		UndoReactionData *data = new UndoReactionData();
		CReactionInterface* ri = new CReactionInterface((*CCopasiRootContainer::getDatamodelList())[0]->getModel());

		if (pModel->getReactions()[i]){
			data->setName(pModel->getReactions()[i]->getObjectName());
			ri->initFromReaction(pModel->getReactions()[i]->getKey());
			data->setRi(ri);
			mpReaData.append(data);
		}
	}
	this->setText(removeAllReactionRowsText());
}

void RemoveAllReactionRowsCommand::redo(){
		mpReactionDM->removeAllReactionRows();
}

void RemoveAllReactionRowsCommand::undo(){
	mpReactionDM->insertReactionRows(mpReaData);
}

QString RemoveAllReactionRowsCommand::removeAllReactionRowsText() const {
	return QObject::tr(": Removed All Reaction Rows");
}

RemoveAllReactionRowsCommand::~RemoveAllReactionRowsCommand() {
	// TODO Auto-generated destructor stub
}

