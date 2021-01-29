/*
    Copyright (c) 2019, Amr Elsawy <aselsawy@eng.asu.edu.eg> or <amrsaad86@gmail.com>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the  Amr Elsawy <amrsaad86@gmail.com> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

        * The use of this code is permitted for research purposes only

    THIS SOFTWARE IS PROVIDED BY Amr Elsawy <amrsaad86@gmail.com> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL Amr Elsawy <amrsaad86@gmail.com> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    Please, cite the following papers if you used this code:

1) Elsawy, Amr S., et al. "MindEdit: A P300-based text editor for mobile devices." Computers in biology and medicine 80 (2017): 97-106.

2) Elsawy, Amr S., et al. "A principal component analysis ensemble classifier for P300 speller applications." 2013 8th International Symposium on Image and Signal Processing and Analysis (ISPA). IEEE, 2013.

3) Elsawy, Amr S., et al. "Performance analysis of a Principal Component Analysis ensemble classifier for Emotiv headset P300 spellers." 2014 36th Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE, 2014.

*/

#ifndef SPELLERGRID_H
#define SPELLERGRID_H

#include <QWidget>
#include <QVector>
#include <QString>
#include "mycallback.h"
#include <hardware/emotiv/sbs2emotivdatareader.h>
#include "overlay.h"


class QLabel;
class QGridLayout;
class QPushButton;
class QVBoxLayout;
class QHBoxLayout;
class QPlainTextEdit;

class SpellerGrid : public QWidget {
    Q_OBJECT

public:

    SpellerGrid (QWidget* parent = 0);

    void closeEvent(QCloseEvent *);
private slots:
    void startThread();
    void stopThread();
    void flashRow(int rowNum);
    void unflashRow(int rowNum);
    void flashCol(int colNum);
    void unflashCol(int colNum);
    void displayTargetCharSlot(QString);
    void startButtonSlot();
    void stopButtonSlot();
    void quitButtonSlot();    
    void finishTrainingSlot();
    void updateOverlay(QString);
    void hideOverlay();
    void emphasizeTarget(QString);
    void demphasizeTarget(QString);

protected:
    void resizeEvent(QResizeEvent *event);

signals:
    void exitSignal();
private:
    int screenHeight;
    int screenWidth;

    static const int nRows = 6;
    static const int nColumns = 6;

    QString Labels;

    QLabel *grid[nRows][nColumns];
    QPlainTextEdit * textEditor;
    QVBoxLayout *vlayout;
    QGridLayout *gridlayout_speller;

    MyCallback* myCallback;
    Sbs2EmotivDataReader* sbs2DataReader;

    //Thread *recordingThread;

    QString targetChar;

    void initEnvironment();
    void initValues();
    void createFolders();
    void createUIs();
    void createLayout();
    void createCharGrid();
    void createConnections();
    void createSequences();

    Overlay * overlay;
};

#endif
