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

/* includes */
#include <QFile>
#include <QStringList>
#include <QTextStream>
#include <QDebug>
#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;

/* paths */
QString PROJECT_FOLDER_TR = "/sdcard/BrainTexEd/train/";
QString FLASHING_FILE_TR = PROJECT_FOLDER_TR+"flashing.txt";
//QString TARGETCHARS_FILE_TR = PROJECT_FOLDER_TR+"targetChar.txt";
QString SIGNALS_FILE_TR = PROJECT_FOLDER_TR+"signals.txt";
QString STIMULUSCODE_FILE_TR = PROJECT_FOLDER_TR+"stimulusCode.txt";

/* definitions */
#define N_CHANNELS 14
#define WINDOW_SIZE 90
#define MovAvg_WIN 13
#define N_ROWCOL 12
#define N_ROWS 6
#define N_COLS 6
#define N_TRIALS 10
#define N_FLASHINGS 12
#define TOTAL_FLASHINGS (N_FLASHINGS*N_TRIALS)
#define THRESHOLD 0.9995

/* data */
// files data
MatrixXd signals_f; // ok
VectorXi flashing_f; // ok
VectorXi stimulusCode_f; // ok

QString charScreen="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"; // ok
QString targetChars="0THE1QUICK2BROWN3FOX4JUMPS5OVER6A7LAZY8DOG9"; // 43 chars
int  nChars=targetChars.size(); // ok


QVector<QVector<MatrixXd> > extractedSignals; // chars x trials x samples x channels // ok
VectorXd labels; // 0s and 1s .. convert to VectorXi if no errors // ok

int nSamples; // ok


MatrixXd singleFeatureVectors;

MatrixXd finalFeatureVectors;

// needed to be saved for online classifier
RowVectorXd dataMean;
MatrixXd pcaEigVectors;
//VectorXd eigValuesMean;
int nPCs;

VectorXd weightVector;

QVector<int> labelsOfTestedData;

/* functions */

void countSamples()
{
    nSamples=0;

    //qDebug() << FLASHING_FILE_TR;

    QFile file(FLASHING_FILE_TR);

    file.open(QFile::ReadOnly);


    if(file.isOpen())
    {
        qDebug() << "flashing file opened ...";
    }

    else
    {
        qDebug() << "flashing file not opened ... ";
    }

    QTextStream stream(&file); // QTextStream removes \n


    while(!stream.readLine().simplified().isNull())
    {
        nSamples++;
    }

    file.close();
}


/* loading */

void loadSignalsData()
{
    qDebug() << "loading signals";

    QFile file(SIGNALS_FILE_TR);

    file.open(QFile::ReadOnly);


    if(file.isOpen())
    {
        qDebug() << "signals file opened ...";
    }

    else
    {
        qDebug() << "signals file not opened ... ";
    }

    QTextStream stream(&file); // QTextStream removes \n

    QString str=stream.readLine().simplified(); // trimmed() will also work

    QStringList list=str.split(' ');


    signals_f.resize(nSamples,list.size()); // or use setZero()

    for(int r=0;r<nSamples;r++)
    {

        list = str.split(' '); // space separator // duplicated for first iteration as I did it above

        for(int c=0;c< list.size();c++)
        {
            signals_f(r,c)=list.at(c).toDouble();
        }

        str = stream.readLine().simplified();

    }

    file.close();

    qDebug() << "loading signals done";
}


void loadFlashingData()
{

    qDebug() << "loading flashing";
    QFile file(FLASHING_FILE_TR);

    file.open(QFile::ReadOnly);


    if(file.isOpen())
    {
        qDebug() << "flashing file opened ...";
    }

    else
    {
        qDebug() << "flashing file not opened ... ";
    }

    QTextStream stream(&file); // QTextStream removes \n

    QString str=stream.readLine().simplified();

    flashing_f.resize(nSamples);

    for(int r=0;r<nSamples;r++)
    {
        flashing_f(r)=str.toInt();
        str=stream.readLine().simplified();
    }

    file.close();

    qDebug() << "loading flashing done";
}


void loadStimulusCodeData()
{
    qDebug() << "loading Stimulus code";

    QFile file(STIMULUSCODE_FILE_TR);

    file.open(QFile::ReadOnly);


    if(file.isOpen())
    {
        qDebug() << "stimulus code file opened ...";
    }

    else
    {
        qDebug() << "stimulus code file not opened ... ";
    }

    QTextStream stream(&file); // QTextStream removes \n

    QString str=stream.readLine().simplified();

    stimulusCode_f.resize(nSamples);

    for(int r=0;r<nSamples;r++)
    {
        stimulusCode_f(r)=str.toInt();
        str=stream.readLine().simplified();
    }

    file.close();
    qDebug() << "loading Stimulus code done";
}


/* preprocessing */

void applyMovingAverageFilter()
{
    qDebug() << "in moving filter";

    // I used window=25 with frequency 240
    // so here with frequency 128 ... I will use about half the window ... 13
    // remember that the window length is odd

    MatrixXd temp;
    RowVectorXd windowAvg;

    int start=MovAvg_WIN/2; // here we will not add 1 as in c++ arrays indexing begins from 0 not 1
    int end=signals_f.rows()-MovAvg_WIN/2; // instead of decrement by 1 ... just remove = from for_loop condition

    qDebug() << signals_f.rows();

    int rowSpan=end-start; // the 1 to be added is already added to end

    qDebug()<< "--- 1 ---";

    qDebug() << rowSpan;
    qDebug() << signals_f.cols();

    temp.setZero(rowSpan,signals_f.cols());

    qDebug()<< "--- 2 ---";

    int row_i=0; // row index

    for(int r=start;r<end;r++)
    {

        windowAvg=signals_f.block(row_i,0,MovAvg_WIN,signals_f.cols()).colwise().mean();

        temp.row(row_i)=windowAvg;

        row_i++;
    }

    qDebug()<< "--- 3 ---";

    // there  is an aliasing problem ... the new values are used in subsequent calculations ...
    // so I made a copy to avoid writing values before finishing ...
    signals_f.block(start,0,rowSpan,signals_f.cols())=temp;

    qDebug() << "out from  moving filter";

}


void applyCARfilter()
{
   qDebug() << "in CAR filter";
    signals_f = signals_f.colwise()- signals_f.rowwise().mean();
   qDebug() << "out from CAR filter";
}

void applyBandPassFilter()
{
    // filter coefficients

    RowVectorXd bn(9);
    bn <<0.000133237244131,0,-0.000532948976523,0,0.000799423464784,0,-0.000532948976523,0,0.000133237244131;

    RowVectorXd an(9);
    an<<0.546374377153003,-4.691514624161844,17.649413866946905,-37.994434212694181,51.190686341642618,-44.200733102520537,23.884541780040749,-7.384334426278459,1.000000000000000;

    // filter

    int N=bn.size();
    int M=an.size();
    int nRows=signals_f.rows();
    int nCols=signals_f.cols();

    // padding zeros to signal
    MatrixXd temp_sig(nRows+N-1,nCols);
    temp_sig.block(0,0,N-1,nCols)=MatrixXd::Zero(N-1,nCols);
    temp_sig.block(N-1,0,nRows,nCols)=signals_f;

    // padding zeros to filtered signal and initialize it to zeros
    MatrixXd filtered(nRows+M-1,nCols);
    filtered.setZero(nRows+M-1,nCols);

    // looping on each sample to calculate its value using filter

    for(int sample=0;sample<nRows;sample++)
    {
        RowVectorXd feed_forward;

        feed_forward = bn*temp_sig.block(sample,0,N,nCols);

        RowVectorXd feed_backward;

        feed_backward = an*filtered.block(sample,0,M,nCols);

        filtered.row(sample+M-1)=feed_forward-feed_backward;
    }

    signals_f=filtered.block(M-1,0,nRows,nCols);
}

MatrixXd zscore(MatrixXd mysig)
{

    int nRows=mysig.rows();
    int nCols=mysig.cols();

    MatrixXd scaled(nRows,nCols);

    for(int c=0;c<nCols;c++)
    {
        VectorXd centered=mysig.col(c).array()-mysig.col(c).mean();
        //double std_dev=sqrt(centered.array().square().sum()/(centered.size()-1));
        double std_dev=centered.norm()/sqrt(centered.size()-1.0);
        scaled.col(c) = centered/std_dev;
    }

    return scaled;
}


/* */

void createData()
{

    qDebug() << "creating data ...";

    // loop here on all data simulataneously and extract based on flashing
    // typically on last 1 (i.e. followed by zero) extract the last 12 samples

    // buffers to hold rows & cols signals in order
    QVector<MatrixXd> rowcol(N_ROWCOL); // 12 represents the row/col flashings

    // initializing buffers with zeros
    for(int rowColNum=0;rowColNum < rowcol.size() ; rowColNum++) // rc rowcol
    {
        rowcol[rowColNum].setZero(WINDOW_SIZE,N_CHANNELS);
    }

    qDebug() << "--- create 1 ---";

    int count=0; // flashings count for 1 char ...
    int totalCount=0; // flashings count for all chars ...
    QVector<int> trialCount(N_ROWCOL,0);

    qDebug() << "--- create 2 ---";

    for(int sample=0; sample < flashing_f.size();sample++)
    {

        int rowcol_i;

        if( ( flashing_f(sample)==0 ) && ( flashing_f(sample+1)==1 )) // catch start of flashings
        {
            // here do extraction
            // I will extract 128 samples(128Hz * 1s) x 14 channels
            count++;
            totalCount++;

            // decrement 1 to scale from 1-12 to 0-11
            // take the next sample as it is the onset of the stimulus
            rowcol_i=stimulusCode_f(sample+1)-1;

            trialCount[rowcol_i]++;

            if(trialCount.at(rowcol_i)<=N_TRIALS)
            {
                // start = i -12 + 1 (as i is counted in the 12)
				// zscore here
                rowcol[rowcol_i]+=zscore(signals_f.block(sample+1,0,WINDOW_SIZE,signals_f.cols())); // catch window starting from next sample
            }
        }


        if(count == TOTAL_FLASHINGS)
        {
            // divide by 15 to get average

            for(int rowColNum=0 ; rowColNum < rowcol.size() ; rowColNum++) // rc rowcol
            {
                rowcol[rowColNum]/=N_TRIALS; // divide by constant = 15 to get average
            }

            // append to extractedSignals after collecting

            extractedSignals.append(rowcol);

            count=0; // reset count

            // reset buffers ...
            for(int rowColNum=0 ; rowColNum < rowcol.size() ; rowColNum++) // rc rowcol
            {
                rowcol[rowColNum].setZero(WINDOW_SIZE,N_CHANNELS);
            }

            for(int rc=0;rc<N_ROWCOL;rc++ )
            {
                trialCount[rc]=0;
            }


        }

        // to exit when all chars are processed ...
        if(totalCount==TOTAL_FLASHINGS*nChars) // trials * #chars (i.e. get #char from file)
        {
            break;
        }

    }

    //signals_f.resize(0,0);
    //flashing_f.resize(0,0);
    //stimulusCode_f.resize(0,0);

    qDebug() << "out from create data ...";
}


void getLabels()
{


    qDebug() << "creating labels ...";
    // initialize labels vector
    labels.setZero(N_FLASHINGS*nChars); // 12 values for each char

    int offset=0; // to point to each 12 vectors start


    for(int i=0 ; i < targetChars.size() ; i++) // loop on chars
    {
        for(int r=0; r < N_ROWS; r++)
        {
            for(int c=0; c < N_COLS ; c++)
            {
                if(targetChars.at(i)==charScreen.at(c+r*N_COLS))
                {
                    // target row
                    labels(offset+r+N_COLS)=1; // +6 to map from 0-5 to 6-11 (indices in c++ array from 0 so 6-11 is actually 7-12)
                    // target col
                    labels(offset+c)=1; // +1 to map from 0-5 to 1-6

                    break;
                }
            }
        }

        offset+=N_ROWCOL;
    }


    qDebug() << "creating labels done ...";
}


void createSingleFeatureVector()
{

    //qDebug() << "single conc. feature vector start";

    // set size to (#Chars*12) rows & (#samples*#channels) cols ...
    singleFeatureVectors.setZero(nChars*N_ROWCOL,WINDOW_SIZE*N_CHANNELS);

    MatrixXd temp;

    int offset = 0; // to point to each 12 vectors start

    for(int char_i=0; char_i < nChars ; char_i++)
    {
        for(int rowColNum=0; rowColNum < extractedSignals.at(char_i).size() ; rowColNum++) // 12 iteration
        {

           temp= extractedSignals.at(char_i).at(rowColNum);

           // resizing matrix into row vector ...
           // size of matrix returns the number of elemnets

           temp.resize(1,WINDOW_SIZE*N_CHANNELS);
           // _OR_USE_ temp.resize(1,extractedSignals.at(char_i).at(rowColNum).size());



           // set the feature veector as row in the matrix of each character
           singleFeatureVectors.row(rowColNum+offset)=temp;

        }

        offset+=N_ROWCOL;
    }

    //qDebug() << "single conc. feature vector end";
}


void decimate(int factor)
{

    // you do not need to use createData before it

    qDebug() << "decimate channel start";

    // set size to (#Chars*12) rows & (#samples*#channels) cols ...
    singleFeatureVectors.setZero(nChars*N_ROWCOL,WINDOW_SIZE*N_CHANNELS);


    int offset = 0; // to point to each 12 vectors start

    for(int char_i=0; char_i < nChars ; char_i++)
    {
        for(int rowColNum=0; rowColNum < extractedSignals.at(char_i).size() ; rowColNum++) // 12 iteration
        {

            MatrixXd temp(WINDOW_SIZE,N_CHANNELS);

            for(int r=0;r<WINDOW_SIZE;r++)
            {
                temp.row(r)= extractedSignals.at(char_i).at(rowColNum).row(r*factor);
            }

           // resizing matrix into row vector ...
           // size of matrix returns the number of elemnets

           temp.resize(1,WINDOW_SIZE*N_CHANNELS);
           // _OR_USE_ temp.resize(1,extractedSignals.at(char_i).at(rowColNum).size());

           // set the feature veector as row in the matrix of each character
           singleFeatureVectors.row(rowColNum+offset)=temp;

        }

        offset+=N_ROWCOL;
    }

    qDebug() << "decimate channel end";
}



void pca(double threshold)
{
    //qDebug() << "PCA start";

    int n = singleFeatureVectors.rows();

    dataMean=singleFeatureVectors.colwise().mean();
    MatrixXd centered = singleFeatureVectors.rowwise() - dataMean;

    MatrixXd cov = ( centered.adjoint() * centered )/(n-1);

    SelfAdjointEigenSolver<MatrixXd> eig(cov);

    VectorXd eigValues = eig.eigenvalues().reverse(); // reverse to make order ascending ...

    int nPCs=1;

    while(( eigValues.head(nPCs).cwiseAbs().sum() / eigValues.cwiseAbs().sum() ) < threshold)
    {
        nPCs++;
    }

    //qDebug() << "nPCs  " << nPCs;

    MatrixXd eigVectors=eig.eigenvectors().rowwise().reverse();
    // centered dos not change so use it as-is ...

    finalFeatureVectors=centered*eigVectors.leftCols(nPCs);

    pcaEigVectors=eigVectors.leftCols(nPCs);

    //qDebug() << "PCA end";
}



void trainFLD()
{
      qDebug() << "train FLD start";

      MatrixXd XP300; // feature vectors of P300 class
      MatrixXd XNP300; // feature vectors of non-P300 class


      // set sizes
      XP300.resize(nChars*2,finalFeatureVectors.cols()); // we have 2 out of 12 P300 signals for each character
      XNP300.resize(nChars*10,finalFeatureVectors.cols()); // we have 10 out of 12 non-P300 signals for each character

      int row_xp300 = 0;
      int row_xnp300 = 0;

      for(int r=0 ; r < finalFeatureVectors.rows() ; r++)
      {
          if(labels(r)==1)
          {
              XP300.row(row_xp300)=finalFeatureVectors.row(r);
              row_xp300++;
          }
          else
          {
              XNP300.row(row_xnp300)=finalFeatureVectors.row(r);
              row_xnp300++;
          }
      }

      // means ... note the means are colwise
      RowVectorXd M1 = XP300.colwise().mean();
      RowVectorXd M2 = XNP300.colwise().mean();

      /** covariance matrices **/
      MatrixXd centered;

      // P300 class
      int N1=XP300.rows();

      centered = XP300.rowwise() - M1;

      MatrixXd Cov1 = ( centered.adjoint() * centered )/(N1-1);

      // non-P300 class
      int N2=XNP300.rows();

      centered = XNP300.rowwise() - M2;

      MatrixXd Cov2 = ( centered.adjoint() * centered )/(N2-1);

      // within class scatter matrix
      MatrixXd Sw=Cov1+Cov2;

      weightVector = Sw.colPivHouseholderQr().solve((M1-M2).transpose());

      qDebug() << "train FLD end";
}


// can be removed as we just want to train

void testFLD()
{
    qDebug() << "test FLD start";

    VectorXd cols;
    VectorXd rows;
    double maxColScore;
    double maxRowScore;
    VectorXd::Index maxRow, maxCol;

    for(int r=0;r<finalFeatureVectors.rows();r+=N_ROWCOL) // step = 12 (12 vectors / char)
    {
        cols = finalFeatureVectors.block(r,0,N_COLS,finalFeatureVectors.cols())*weightVector;
        maxColScore =cols.maxCoeff(&maxCol);

        rows = finalFeatureVectors.block(r+N_COLS,0,N_ROWS,finalFeatureVectors.cols())*weightVector;
        maxRowScore=rows.maxCoeff(&maxRow);

        labelsOfTestedData.append(N_COLS*maxRow+maxCol);
    }

    qDebug() << "test FLD end";
}


// can be removed as we just want to train
void printText()
{
    qDebug() << "print text start";

//    QString str="A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,0,1,2,3,4,5,6,7,8,9";
//    QStringList charScreen = str.split(',');

    QTextStream out(stdout);
    for(int i=0;i<labelsOfTestedData.size();i++)
    {
        out << charScreen.at(labelsOfTestedData.at(i));
    }

    qDebug() << "print text end";
}


void saveWeightVector()
{
    QFile output("/sdcard/BrainTexEd/par_single/weightVector");
    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    // save size at first line
    stream << weightVector.size() << endl;

    for(int i=0;i<weightVector.size();i++)
    {
            stream << weightVector(i) << endl;
    }

    output.close();

}


void saveDataMean()
{
    QFile output("/sdcard/BrainTexEd/par_single/dataMean");
    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    // save size at first line
    stream << dataMean.size() << endl;

    for(int i=0;i<dataMean.size();i++)
    {
            stream << dataMean(i) << endl;
    }

    output.close();

}


void saveEigVectors()
{
    QFile output("/sdcard/BrainTexEd/par_single/eigVectors");
    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    // save size at first two lines
    stream << pcaEigVectors.rows() << endl;
    stream << pcaEigVectors.cols() << endl;

    for(int r=0;r<pcaEigVectors.rows();r++)
    {
        for(int c=0;c<pcaEigVectors.cols();c++)
        {
            stream << pcaEigVectors(r,c) << " ";
        }

        stream << endl;
    }

    output.close();
}


void train_eigen_single()
{
    // list all operations in order
    countSamples();
    getLabels();
    loadSignalsData();
    applyCARfilter();
    //applyMovingAverageFilter();
	applyBandPassFilter();
    loadFlashingData();
    loadStimulusCodeData();
    createData();
    createSingleFeatureVector();
    pca(THRESHOLD);
    trainFLD();
    testFLD();
    printText();
    saveWeightVector();
    saveDataMean();
    saveEigVectors();

}
