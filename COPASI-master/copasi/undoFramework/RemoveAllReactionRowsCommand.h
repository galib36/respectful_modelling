/*
 * RemoveAllReactionRowsCommand.h
 *
 *  Created on: 12 Aug 2014
 *      Author: dada
 */

#ifndef REMOVEALLREACTIONROWSCOMMAND_H_
#define REMOVEALLREACTIONROWSCOMMAND_H_

#include "CCopasiUndoCommand.h"

class RemoveAllReactionRowsCommand: public CCopasiUndoCommand {
public:
	RemoveAllReactionRowsCommand(CQReactionDM * pReaDM, const QModelIndex&);
	void redo();
	void undo();
	QString removeAllReactionRowsText() const;
	virtual ~RemoveAllReactionRowsCommand();

private:
	CQReactionDM *mpReactionDM;
	QList <UndoReactionData *> mpReaData;
};

#endif /* REMOVEALLREACTIONROWSCOMMAND_H_ */
