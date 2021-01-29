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


#include "mycallback.h"
#include <QFile>
#include <QTextStream>
#include <QDir>
#include <QMutex>

#define nRowCols 12
#define nRows 6
#define nCols 6
#define nTrials 10

MyCallback::MyCallback(QObject *parent) :
    Sbs2Callback(parent)
{

    flashingCount=0; //13 -> 100 ms
    unflashingCount=0; //29 -> 225 ms
    beforeEmphasize=0; //128 -> 1 sec
    afterEmphasize=0; //192 -> 1.5 sec // 39 -> 300 ms
    charCounter=0;
    randomSeqCounter=0;
    finishCounter=0;
    totalRepeat=0; // 16 sym * 10 trial
    breaksCount=0;
    totalBreak=1;

    mySlot=MyCallback::preEmphasize;

    recordingFlag=0;
    flashingIndex=0;
    stimulusCode=0;


    // creating file for signals data ...
    signalsFile = new QFile("/sdcard/BrainTexEd/train/signals.txt");
    signalsFile->open(QFile::WriteOnly);

    // creating file for falshing data ...
    flashingFile = new QFile("/sdcard/BrainTexEd/train/flashing.txt");
    flashingFile->open(QFile::WriteOnly);


    // creating file for stimulusCode data ...
    stimulusCodeFile = new QFile("/sdcard/BrainTexEd/train/stimulusCode.txt");
    stimulusCodeFile->open(QFile::WriteOnly);


    QString str="T,H,E,Q,U,I,C,K,B,R,O,W,N,F,O,X,J,U,M,P,S,O,V,E,R,A,L,A,Z,Y,D,O,G";

    chars=str.split(',');

    nChars=chars.size();


    printChannels();

    createRandomSequence();
    printRandSequences();

}



// getData() is called from execute() method in sbs2eomtivdatareader in hardware/emotiv
// it is called in ifinite loop
// each time it takes one packet ...
// so this place is the potential place for recording other speller data (ex flashing ... stimulus code ... )
// make sure that packet represent one readings from all channels
// then for each packet make
// flag: flashing = 1 or 0
// stimulsCode = 0,1,...,12
//
void MyCallback::getData(Sbs2Packet *packet)
{
    // passing packet to callback and datahandler objects ...
    setPacket(packet);

    // applying filter ... it will execute based on whether the filter is turned on or not
    //sbs2DataHandler->filter();


    if(recordingFlag)
    {


       if(mySlot==breaks)
       {
           breaksCount++;

           if(breaksCount==totalBreak)
           {
                breaksCount=0;
                mySlot=MyCallback::preEmphasize;
           }
       }

       else if(mySlot==preEmphasize)
       {

           beforeEmphasize++;

          if(beforeEmphasize==128)
          {
              beforeEmphasize=0;
              mySlot=MyCallback::postEmphasize;
              // should be emitted each app flashing
              emit breakEndSignal();
              emit emphasizeTargetSignal(chars.at(charCounter));
              emit displayTargetSignal(chars.at(charCounter));
          }

          recordall(packet);
       }

       else if(mySlot==MyCallback::postEmphasize)
       {

           afterEmphasize++;

            if(afterEmphasize==192)
            {
                afterEmphasize=0;

                mySlot=MyCallback::flashing;

                emit demphasizeTargetSignal(chars.at(charCounter));

            }

            recordall(packet);
       }

       else if(mySlot==MyCallback::flashing)
       {

           int rowORcol=randSequences.at(randomSeqCounter);

           if(flashingCount==0)
           {

              if(rowORcol>nCols)// nCols=6 where rows begin from 7
              {
                  // signal that makes flashing begins
                  emit rowFlashBeginSignal(rowORcol);
              }
              else
              {
                  // signal that makes flashing begins
                  emit colFlashBeginSignal(rowORcol);
              }
           }

          flashingIndex=1;
          stimulusCode=rowORcol;

            flashingCount++;

            if(flashingCount==13)
            {
                flashingCount=0;

                mySlot=MyCallback::unflashing;
            }

            recordall(packet);
       }

       else if(mySlot==MyCallback::unflashing)
       {

           int rowORcol=randSequences.at(randomSeqCounter);

          if (unflashingCount==0)
          {
              if(rowORcol>nCols)// nCols=6 where rows begin from 7
              {
                  // signal that makes flashing begins
                  emit rowFlashEndSignal(rowORcol);
              }
              else
              {
                  // signal that makes flashing begins
                  emit colFlashEndSignal(rowORcol);
              }
          }

           flashingIndex=0;
           stimulusCode=0;

           unflashingCount++;

            if(unflashingCount==39)// 10 29 39
            {
                unflashingCount=0;
                mySlot=MyCallback::flashing;
                randomSeqCounter++;
                totalRepeat++;

                if(totalRepeat==nRowCols*nTrials)
                {
                    totalRepeat=0;
                    mySlot=MyCallback::breaks;

                    charCounter++;

                    if(charCounter==10||charCounter==20||charCounter==30||charCounter==40)
                    {
                        totalBreak=60*128; // 1 min break
                        emit breakBeginSignal("break");
                    }
                    else
                    {
                        totalBreak=1;
                    }

                    if (charCounter==nChars)
                    {
                         mySlot=MyCallback::finish;
                    }
                }

            }

            recordall(packet);
       }

       else if(mySlot==MyCallback::finish)
       {
           finishCounter++;

           if(finishCounter==128)
           {
               recordingFlag=0;
           }

           recordall(packet);
       }


       // record


       // do training here
       if(recordingFlag==0)
       {

            emit finishTrainingSignal(); //

       }

    }


    if (!(thisPacket->cq == -1))
       emit cqValues(thisPacket->cqName,thisPacket->cq);


    int packetInRecording = (currentPacket - sbs2DataHandler->getPacketZero() +1);
    if (packetInRecording%4 == 0)
    emit currentPacketInRecording(packetInRecording);

}



void MyCallback::recordall(Sbs2Packet *packet)
{
    // recording into signals txt file ...
    recordSignals(packet);

    // record flashing ...
    recordFlashing();
    recordStimulusCode();
}

void MyCallback::setRecordingFlag()
{
    recordingFlag=1;
}

void MyCallback::resetRecordingFlag()
{
    recordingFlag=0;
}

void MyCallback::turnOnFilter()
{
    turnFilterOn(8,30,32);
}


void MyCallback::recordFlashing()
{

    QTextStream toFlashingFile(flashingFile);
    toFlashingFile << flashingIndex << endl ;
}

void MyCallback::recordSignals(Sbs2Packet *packet)
{
    // I did not use the reference & as I am using a pointer
    QTextStream toSignalsFile(signalsFile);


    // save these values to file in one line
    for (int i=0; i<Sbs2Common::channelsNo(); ++i)
    {
        toSignalsFile << ( packet->filteredValues[Sbs2Common::getChannelNames()->at(i)] ) << " ";
    }

        toSignalsFile << endl ;
}

void MyCallback::printChannels()
{

    for (int i=0; i<Sbs2Common::channelsNo(); ++i)
    {
        qDebug() << Sbs2Common::getChannelNames()->at(i);
    }

}

void MyCallback::recordStimulusCode()
{
    QTextStream toStimulusCodeFile(stimulusCodeFile);
    toStimulusCodeFile << stimulusCode << endl ;
}



MyCallback::~MyCallback()
{
    // close signals file ...
    signalsFile->close();
    flashingFile->close();
    stimulusCodeFile->close();

//    QString _signals = QString("/sdcard/BrowserData/br_bpf_online/signals-%1.txt").arg(timestamp);
//    QString _flashing = QString("/sdcard/BrowserData/br_bpf_online/flashing-%1.txt").arg(timestamp);
//    QString _stimulusCode = QString("/sdcard/BrowserData/br_bpf_online/stimulusCode-%1.txt").arg(timestamp);

//    QFile::rename("/sdcard/BrowserData/br_bpf_online/signals.txt",_signals);
//    QFile::rename("/sdcard/BrowserData/br_bpf_online/flashing.txt",_flashing);
//    QFile::rename("/sdcard/BrowserData/br_bpf_online/stimulusCode.txt",_stimulusCode);

}



void MyCallback::createRandomSequence()
{
    QVector<int> seedSequence;

    int N=nRowCols;

    // initialize seed sequence with numbers from 1 to 12
    for(int i=1 ; i<= N ; i++)
    {
        seedSequence.append(i);
    }

    int pos,randNum,temp;

    for(int trial=0 ; trial < nTrials*nChars ; trial++)
    {

        // shuffle

         pos=N-1; // points to last element in seed Sequence

        for(int i=0 ; i< N-1 ; i++)
        {
            randNum=rand()%(pos);

            temp=seedSequence[pos];
            seedSequence[pos]=seedSequence[randNum];
            seedSequence[randNum]=temp;

            pos--;

        }

        // transfer shuffled numbers to the whole randSequences vector

        for(int i=0 ; i< N ; i++)
        {
            randSequences.append(seedSequence.at(i));
        }

    }

}

void MyCallback::printRandSequences()
{

    int N=nRowCols;

    for(int i=0 ; i< randSequences.size() ; i++)
    {
        if(i%N==0) std::cout << std::endl;
        std::cout << randSequences.at(i) << " ";

    }
}

