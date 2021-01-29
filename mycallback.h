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

#ifndef MYCALLBACK_H
#define MYCALLBACK_H

#include <sbs2callback.h>
#include <sbs2datahandler.h>
class QFile;


class MyCallback : public Sbs2Callback
{
    Q_OBJECT
public:
    explicit MyCallback(QObject *parent = 0);//,int session
     ~MyCallback();
    void getData(Sbs2Packet *packet);

signals:
    void currentPacketInRecording(QVariant value);
    void displayTargetSignal(QString);
    void emphasizeTargetSignal(QString);
    void demphasizeTargetSignal(QString);
    void rowFlashBeginSignal(int);
    void rowFlashEndSignal(int);
    void colFlashBeginSignal(int);
    void colFlashEndSignal(int);
    void finishTrainingSignal();
    void breakBeginSignal(QString);
    void breakEndSignal();

public slots:
    void turnOnFilter();
    void setRecordingFlag();
    void resetRecordingFlag();



private:
    QFile *signalsFile;
    QFile *flashingFile;
    QFile *stimulusCodeFile;
    QFile *targetCharFile;


    int flashingCount;
    int unflashingCount;
    int beforeEmphasize;
    int afterEmphasize;
    int freeCounter;
    int charCounter;
    int randomSeqCounter;
    int finishCounter;
    int totalRepeat;
    int breaksCount;
    int totalBreak;

    int recordingFlag;

    int flashingIndex;
    int stimulusCode;

    void recordFlashing();
    void recordSignals(Sbs2Packet *packet);
    void recordStimulusCode();
    void recordTargetChar();
    void printChannels();
    void recordall(Sbs2Packet *packet);

    enum TimeSlot{
        breaks,preEmphasize,postEmphasize,flashing,unflashing,finish
    } mySlot;

    QVector<int> randSequences;

    QStringList chars;
    int nChars;

    void createRandomSequence();
    void printRandSequences();
};

#endif // MYCALLBACK_H
