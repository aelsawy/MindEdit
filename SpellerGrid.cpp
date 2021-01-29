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

    THIS SOFTWARE IS PROVIDED BY Amr Elsawt <amrsaad86@gmail.com> ''AS IS'' AND ANY
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

#include "SpellerGrid.h"
#include "cstdlib"
#include <QLabel>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QPlainTextEdit>
#include <QDebug>
#include <QString>
#include <QPalette>
#include <QColor>
#include <QTextEdit>
#include <QFile>
#include <QTextStream>

#include "trainer.h"
#include "trainer_linear.h"


SpellerGrid::SpellerGrid (QWidget* parent) :
    QWidget (parent)
    {
        createFolders();
        initEnvironment();
        initValues();
        createUIs();
        createLayout();
        createConnections();
        startButtonSlot();

        overlay = new Overlay(this);
        overlay->hide();
        //overlay->setShown(false);

    }


void SpellerGrid::resizeEvent(QResizeEvent *event) {
    // size using ratio of widget size
    overlay->resize(this->width()/2,this->height()/2);
    overlay->move(this->width()/4,this->height()/4);
    event->accept();
}

void SpellerGrid::initEnvironment()
{
    qDebug() << "catalogPath: "<<Sbs2Common::setDefaultCatalogPath();
    qDebug() << "rootAppPath: "<<Sbs2Common::setDefaultRootAppPath();

    // I added this pointer so check for whether it is ok
    myCallback = new MyCallback(this);
    sbs2DataReader = Sbs2EmotivDataReader::New(myCallback,0);

    Sbs2Common::setHardware("emotiv");

}


void SpellerGrid::initValues()
{

    Labels="A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,0,1,2,3,4,5,6,7,8,9";

    //recordingThread = new Thread(this);
    //recordingThread->session=session;
    //this->showMaximized();
    //this->showFullScreen();


}


void SpellerGrid::closeEvent(QCloseEvent *event)
{
    this->quitButtonSlot();
    //event->accept();
}



void SpellerGrid::createFolders()
{
    QDir().mkdir("/sdcard/BrainTexEd");
    QDir().setCurrent("/sdcard/BrainTexEd");
    QDir().mkdir("/sdcard/BrainTexEd/par_ensemble");
    QDir().mkdir("/sdcard/BrainTexEd/par_linear");
    QDir().mkdir("/sdcard/BrainTexEd/par_single");
    QDir().mkdir("/sdcard/BrainTexEd/train");

}


void SpellerGrid::createConnections()
{

    QObject::connect(myCallback,SIGNAL(displayTargetSignal(QString)),this,SLOT(displayTargetCharSlot(QString)));

    QObject::connect(myCallback, SIGNAL(rowFlashBeginSignal(int)), this, SLOT(flashRow(int)));
    QObject::connect(myCallback, SIGNAL(colFlashBeginSignal(int)), this, SLOT(flashCol(int)));

    QObject::connect(myCallback, SIGNAL(rowFlashEndSignal(int)), this, SLOT(unflashRow(int)));
    QObject::connect(myCallback, SIGNAL(colFlashEndSignal(int)), this, SLOT(unflashCol(int)));

    QObject::connect(myCallback, SIGNAL(emphasizeTargetSignal(QString)), this, SLOT(emphasizeTarget(QString)));
    QObject::connect(myCallback, SIGNAL(demphasizeTargetSignal(QString)), this, SLOT(demphasizeTarget(QString)));



    QObject::connect(myCallback, SIGNAL(breakBeginSignal(QString)), this, SLOT(updateOverlay(QString)));
    QObject::connect(myCallback, SIGNAL(breakEndSignal()), this, SLOT(hideOverlay()));

    QObject::connect(myCallback, SIGNAL(finishTrainingSignal()), this, SLOT(finishTrainingSlot()));


    // speller exit symbol signal
    QObject::connect(this, SIGNAL(exitSignal()), this, SLOT(quitButtonSlot()));

    QObject::connect(textEditor, SIGNAL(cursorPositionChanged()), this, SLOT(highlightLineSlot()));
}




void SpellerGrid::emphasizeTarget(QString targetChar)
{

    qDebug()<<"emphasizing entry ...";


    QString str="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

    int index=str.indexOf(targetChar);

    int r=index/nColumns;
    int c=index%nColumns;

    grid[r][c]->setStyleSheet("QLabel { color : red; }");
    QFont cellFont( "Arial", 40 , QFont::Bold);
    grid[r][c]->setFont( cellFont);

    qDebug()<<"emphasizing exit ...";

}


void SpellerGrid::demphasizeTarget(QString targetChar)
{

    qDebug()<<"demphasizing entry ...";

    QString str="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

    int index=str.indexOf(targetChar);

    int r=index/nColumns;
    int c=index%nColumns;

    grid[r][c]->setStyleSheet("QLabel { color : #222; }");
    QFont cellFont( "Arial", 30 , QFont::Bold);
    grid[r][c]->setFont( cellFont);

    qDebug()<<"demphasizing exit ...";

}

void SpellerGrid::finishTrainingSlot()
{

    QString path="/sdcard/BrainTexEd/Time.txt";
    QFile output(path);

    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    QElapsedTimer timer;
    timer.start();

    //
    train_eigen();

    stream << timer.restart() << endl;

    //
    train_linear();

    stream << timer.restart() << endl;

    //
    train_eigen_single();

    stream << timer.elapsed() << endl;

    output.close();

    updateOverlay("End");

}


void SpellerGrid::updateOverlay(QString msg)
{
    overlay->setText(msg);
    overlay->repaint();
    overlay->show();
}


void SpellerGrid::hideOverlay()
{
    overlay->hide();
}

void SpellerGrid::startButtonSlot()
{
    qDebug()<<"start button slot entry ...";

    QMutex mutex;
    mutex.lock();

    myCallback->setRecordingFlag();
    myCallback->turnOnFilter();

    // I need to clear text ...
    // close old files
    // open new files

    //this->startThread();

    textEditor->setFocus();

    mutex.unlock();

    qDebug()<<"start button slot exit ...";
}

void SpellerGrid::stopButtonSlot()
{
    qDebug()<<"stop button slot entry ...";

    QMutex mutex;
    mutex.lock();

    myCallback->resetRecordingFlag();
    myCallback->turnFilterOff();
    //this->stopThread();

    mutex.unlock();

    qDebug()<<"stop button slot exit ...";
}

void SpellerGrid::quitButtonSlot()
{
    qDebug()<<"quit button slot entry ...";

    QMutex mutex;
    mutex.lock();

    myCallback->resetRecordingFlag();
    myCallback->turnFilterOff();
      //this->stopThread();
    sbs2DataReader->aboutToQuit();
      this->close();

    mutex.unlock();

    qDebug()<<"quit button slot exit ...";
}


void SpellerGrid::displayTargetCharSlot(QString targetChar)
{

    QMutex mutex;
    mutex.lock();

    qDebug()<<"displaying target entry ...";

    this->targetChar=targetChar;
    textEditor->insertPlainText(targetChar);
    textEditor->ensureCursorVisible();

    qDebug()<<"displaying target exit ...";

    mutex.unlock();
}


void SpellerGrid::startThread()
{
    //recordingThread->stop=false;
    //recordingThread->start();
}

void SpellerGrid::stopThread()
{
    //recordingThread->stop=true;
}


void SpellerGrid::createUIs()
{
    textEditor = new QPlainTextEdit();
    QFont textFont( "Arial",20);//15
    textEditor->setFont(textFont);
    textEditor->setStyleSheet("QPlainTextEdit { background-color : #222; color : white;}");
    //textEditor->setCursorWidth(10);
    textEditor->setCenterOnScroll(true);
    textEditor->setEnabled(false);
    textEditor->setFixedHeight(200);
    textEditor->setWordWrapMode();

//    QTextEdit::ExtraSelection highlight;
//    highlight.cursor = textEditor->textCursor();
//    highlight.format.setProperty(QTextFormat::FullWidthSelection, false);
//    highlight.format.setBackground(QColor(202,255,255));//QColor(30,144,255)
//    QList<QTextEdit::ExtraSelection> extras;
//    extras << highlight;
//    textEditor->setExtraSelections( extras );

    createCharGrid();
}


void SpellerGrid::createLayout()
{

    vlayout = new QVBoxLayout;
    vlayout->addWidget(textEditor);
    vlayout->addLayout(gridlayout_speller);

    this->setLayout(vlayout);
}

void SpellerGrid::createCharGrid()
{
    QString screenLabels="A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,0,1,2,3,4,5,6,7,8,9";
    QStringList list=screenLabels.split(',');

    gridlayout_speller= new QGridLayout;

    for (int r = 0; r < nRows; r++) {

        for (int c = 0; c < nColumns; c++) {

            QLabel* cell = new QLabel(list.at(c+nColumns*r));

            // cell font type and size
            QFont cellFont( "Arial", 25 , QFont::Bold);//35
            cell->setFont( cellFont);

            // each cell text color and background color
            cell->setStyleSheet("QLabel { background-color : black; color : #222; }");
            //cell->setFixedHeight(165);
            cell->setFixedSize(120,165);
            //
            //cell->setStyleSheet("QLabel { qproperty-alignment: AlignCenter; }");
            cell->setAlignment(Qt::AlignCenter);
            cell->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
            // add to grid
            grid[r][c] = cell;
            gridlayout_speller->addWidget(cell, r, c);


        }
    }    

    // QWidget Backrgound color
    setStyleSheet("background-color:black;");
}



void SpellerGrid::flashRow(int row)
{

    QMutex mutex;
    mutex.lock();

    int r=row;

    r=r-(nColumns+1); // we sub 7 as we map from 7-12 to 0-5 (i.e. array indices)

    for (int c = 0; c < nColumns; c++)
    {
        // each cell text color and background color
        grid[r][c]->setStyleSheet("QLabel { color : white; }");
        QFont cellFont( "Arial", 32 , QFont::Bold);
        grid[r][c]->setFont( cellFont);

    }

    mutex.unlock();

    qDebug() << "flash row done";
}


void SpellerGrid::unflashRow(int r)
{

    QMutex mutex;
    mutex.lock();

    r=r-(nColumns+1);

    for (int c = 0; c < nColumns; c++)
    {
        // each cell text color and background color
        grid[r][c]->setStyleSheet("QLabel { color : #222; }");
        QFont cellFont( "Arial", 25 , QFont::Bold);//35
        grid[r][c]->setFont( cellFont);

    }

    mutex.unlock();

    qDebug() << "unflash row done";
}


void SpellerGrid::flashCol(int col)
{

    QMutex mutex;
    mutex.lock();

    int c= col;

    c--; // to map from 1-6 to 0-5

    for (int r = 0; r < nRows; r++)
    {
        // each cell text color and background color
        grid[r][c]->setStyleSheet("QLabel { color : white; }");
        QFont cellFont( "Arial", 32 , QFont::Bold);
        grid[r][c]->setFont( cellFont);

    }

    mutex.unlock();

    qDebug() << "flash col done";
}


void SpellerGrid::unflashCol(int c)
{

    QMutex mutex;
    mutex.lock();

    c--;

    for (int r = 0; r < nRows; r++)
    {
        // each cell text color and background color
        grid[r][c]->setStyleSheet("QLabel { color :#222; }");
        QFont cellFont( "Arial", 25 , QFont::Bold);//35
        grid[r][c]->setFont( cellFont);

    }


    mutex.unlock();

    qDebug() << "unflash col done";
}
