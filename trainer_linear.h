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
QString PROJECT_FOLDER = "/sdcard/BrainTexEd/train/";
QString FLASHING_FILE = PROJECT_FOLDER+"flashing.txt";
//QString TARGETCHARS_FILE = PROJECT_FOLDER+"targetChar.txt";
QString SIGNALS_FILE = PROJECT_FOLDER+"signals.txt";
QString STIMULUSCODE_FILE = PROJECT_FOLDER+"stimulusCode.txt";

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
MatrixXd signalz_m; // ok
VectorXi flashing_v; // ok
VectorXi stimulusCode_v; // ok

QString charzScreen="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"; // ok
QString targetCharz="THEQUICKBROWNFOXJUMPSOVERALAZYDOG"; // 33 chars //"AAAAAAAAAA";
int  nCharz=targetCharz.size(); // ok


QVector<QVector<MatrixXd> > extractedSignalz; // chars x trials x samples x channels // ok
VectorXd labelz; // 0s and 1s .. convert to VectorXi if no errors // ok

int nSamplez; // ok


QVector<MatrixXd> featureVectorzEnsemble; // ok

QVector<MatrixXd> featureVectorz; // by applying PCA or PCA+decimation // ok // featureMatrix // featureVectorz // reducedFeatureVectorz

// needed to be saved for online classifier
QVector<RowVectorXd> ensembleMeanz; // ok
QVector<MatrixXd> ensembleEigVectorz; // ok
VectorXd eigValuezMean; // ok
int nPCz; // ok

QVector<VectorXd> weightVectorz; // to hold all weight vectors // ok

QVector<int> labelzOfTestedData; // ok


// level two classifier
MatrixXd trainingScores;
VectorXd Weight; // weight
//double bias; // bias


/* functions */

void countSamplez()
{
    nSamplez=0;

    //qDebug() << FLASHING_FILE;

    QFile file(FLASHING_FILE);

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
        nSamplez++;
    }

    file.close();
}


/* loading */

void loadSignalzData()
{
    qDebug() << "loading signalz";

    QFile file(SIGNALS_FILE);

    file.open(QFile::ReadOnly);


    if(file.isOpen())
    {
        qDebug() << "signalz file opened ...";
    }

    else
    {
        qDebug() << "signalz file not opened ... ";
    }

    QTextStream stream(&file); // QTextStream removes \n

    QString str=stream.readLine().simplified(); // trimmed() will also work

    QStringList list=str.split(' ');


    signalz_m.resize(nSamplez,list.size()); // or use setZero()

    for(int r=0;r<nSamplez;r++)
    {

        list = str.split(' '); // space separator // duplicated for first iteration as I did it above

        for(int c=0;c< list.size();c++)
        {
            signalz_m(r,c)=list.at(c).toDouble();
        }

        str = stream.readLine().simplified();

    }

    file.close();

    qDebug() << "loading signals done";
}


void loadFlashingVector()
{

    qDebug() << "loading flashing vector";
    QFile file(FLASHING_FILE);

    file.open(QFile::ReadOnly);


    if(file.isOpen())
    {
        qDebug() << "flashing vector file opened ...";
    }

    else
    {
        qDebug() << "flashing vector file not opened ... ";
    }

    QTextStream stream(&file); // QTextStream removes \n

    QString str=stream.readLine().simplified();

    flashing_v.resize(nSamplez);

    for(int r=0;r<nSamplez;r++)
    {
        flashing_v(r)=str.toInt();
        str=stream.readLine().simplified();
    }

    file.close();

    qDebug() << "loading flashing vector done";
}


void loadStimulusCodeVector()
{
    qDebug() << "loading Stimulus code vector";

    QFile file(STIMULUSCODE_FILE);

    file.open(QFile::ReadOnly);


    if(file.isOpen())
    {
        qDebug() << "stimulus code vector file opened ...";
    }

    else
    {
        qDebug() << "stimulus code vector file not opened ... ";
    }

    QTextStream stream(&file); // QTextStream removes \n

    QString str=stream.readLine().simplified();

    stimulusCode_v.resize(nSamplez);

    for(int r=0;r<nSamplez;r++)
    {
        stimulusCode_v(r)=str.toInt();
        str=stream.readLine().simplified();
    }

    file.close();
    qDebug() << "loading Stimulus code vector done";
}


/* preprocessing */

void applyMovingAverageFilterz()
{
    qDebug() << "in moving filterz";

    // I used window=25 with frequency 240
    // so here with frequency 128 ... I will use about half the window ... 13
    // remember that the window length is odd

    MatrixXd temp;
    RowVectorXd windowAvg;

    int start=MovAvg_WIN/2; // here we will not add 1 as in c++ arrays indexing begins from 0 not 1
    int end=signalz_m.rows()-MovAvg_WIN/2; // instead of decrement by 1 ... just remove = from for_loop condition

    qDebug() << signalz_m.rows();

    int rowSpan=end-start; // the 1 to be added is already added to end

    qDebug()<< "--- 1 ---";

    qDebug() << rowSpan;
    qDebug() << signalz_m.cols();

    temp.setZero(rowSpan,signalz_m.cols());

    qDebug()<< "--- 2 ---";

    int row_i=0; // row index

    for(int r=start;r<end;r++)
    {

        windowAvg=signalz_m.block(row_i,0,MovAvg_WIN,signalz_m.cols()).colwise().mean();

        temp.row(row_i)=windowAvg;

        row_i++;
    }

    qDebug()<< "--- 3 ---";

    // there  is an aliasing problem ... the new values are used in subsequent calculations ...
    // so I made a copy to avoid writing values before finishing ...
    signalz_m.block(start,0,rowSpan,signalz_m.cols())=temp;

    qDebug() << "out from  moving filterz";

}


void applyCARfilterz()
{
   qDebug() << "in CAR filterz";
    signalz_m = signalz_m.colwise()- signalz_m.rowwise().mean();
   qDebug() << "out from CAR filterz";
}

void applyBandPassFilterz()
{


    qDebug() << "in BPF filterz";

    // filter coefficients

    RowVectorXd bn(9);
    bn <<0.000133237244131,0,-0.000532948976523,0,0.000799423464784,0,-0.000532948976523,0,0.000133237244131;

    RowVectorXd an(9);
    an<<0.546374377153003,-4.691514624161844,17.649413866946905,-37.994434212694181,51.190686341642618,-44.200733102520537,23.884541780040749,-7.384334426278459,1.000000000000000;

    // filter

    int N=bn.size();
    int M=an.size();
    int nRows=signalz_m.rows();
    int nCols=signalz_m.cols();

    // padding zeros to signal
    MatrixXd temp_sig(nRows+N-1,nCols);
    temp_sig.block(0,0,N-1,nCols)=MatrixXd::Zero(N-1,nCols);
    temp_sig.block(N-1,0,nRows,nCols)=signalz_m;

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

    signalz_m=filtered.block(M-1,0,nRows,nCols);

    qDebug() << "out of BPF filterz";
}

MatrixXd zscorez(MatrixXd mysig)
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

void createDataz()
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

    for(int sample=0; sample < flashing_v.size();sample++)
    {

        int rowcol_i;

        if( ( flashing_v(sample)==0 ) && ( flashing_v(sample+1)==1 )) // catch start of flashings
        {
            // here do extraction
            // I will extract 128 samples(128Hz * 1s) x 14 channels
            count++;
            totalCount++;

            // decrement 1 to scale from 1-12 to 0-11
            // take the next sample as it is the onset of the stimulus
            rowcol_i=stimulusCode_v(sample+1)-1;

            trialCount[rowcol_i]++;

            if(trialCount.at(rowcol_i)<=N_TRIALS)
            {
                // start = i -12 + 1 (as i is counted in the 12)
				// zscore here
                rowcol[rowcol_i]+=zscorez(signalz_m.block(sample+1,0,WINDOW_SIZE,signalz_m.cols())); // catch window starting from next sample
            }
        }


        if(count == TOTAL_FLASHINGS)
        {
            // divide by 15 to get average

            for(int rowColNum=0 ; rowColNum < rowcol.size() ; rowColNum++) // rc rowcol
            {
                rowcol[rowColNum]/=N_TRIALS; // divide by constant = 15 to get average
            }

            // append to extractedSignalz after collecting

            extractedSignalz.append(rowcol);

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
        if(totalCount==TOTAL_FLASHINGS*nCharz) // trials * #chars (i.e. get #char from file)
        {
            break;
        }

    }

    //signalz_m.resize(0,0);
    //flashing_v.resize(0,0);
    //stimulusCode_v.resize(0,0);

    qDebug() << "out from create data ...";
}


void getLabelz()
{


    qDebug() << "creating labelz ...";
    // initialize labelz vector
    labelz.setZero(N_FLASHINGS*nCharz); // 12 values for each char

    int offset=0; // to point to each 12 vectors start


    for(int i=0 ; i < targetCharz.size() ; i++) // loop on chars
    {
        for(int r=0; r < N_ROWS; r++)
        {
            for(int c=0; c < N_COLS ; c++)
            {
                if(targetCharz.at(i)==charzScreen.at(c+r*N_COLS))
                {
                    // target row
                    labelz(offset+r+N_COLS)=1; // +6 to map from 0-5 to 6-11 (indices in c++ array from 0 so 6-11 is actually 7-12)
                    // target col
                    labelz(offset+c)=1; // +1 to map from 0-5 to 1-6

                    break;
                }
            }
        }

        offset+=N_ROWCOL;
    }


    qDebug() << "creating labelz done ...";
}


void createEnsembleFeatureVectorz()
{

    qDebug() << "single Ensemble feature vector start";

    qDebug() << "samples = " << WINDOW_SIZE << " channels = " << N_CHANNELS;
    featureVectorzEnsemble.resize(N_CHANNELS); // 14 matrices


    // set size to (#Chars*12) rows & (#samples) cols ... for each matrix
    for(int channel=0 ; channel < N_CHANNELS ; channel++)
    {
        featureVectorzEnsemble[channel].setZero(nCharz*N_ROWCOL,WINDOW_SIZE);
    }

    qDebug() << "r&c " << featureVectorzEnsemble.at(0).rows() << " " << featureVectorzEnsemble.at(0).cols();

    MatrixXd temp;

    int offset = 0; // to point to each 12 vectors start

    for(int char_i=0; char_i < nCharz ; char_i++)
    {
        for(int rowColNum=0; rowColNum < extractedSignalz.at(char_i).size() ; rowColNum++) // 12 iteration
        {

           temp= extractedSignalz.at(char_i).at(rowColNum);

           for(int channel=0 ; channel < N_CHANNELS ; channel++)
           {

                // set the col as feature vector row in the matrix of each channel of each character
               featureVectorzEnsemble[channel].row(rowColNum+offset)=temp.col(channel).transpose();
           }


        }

        offset+=N_ROWCOL;
    }

    qDebug() << "single Ensemble feature vector end";
}

void pcaz(double threshold)
{
    qDebug() << "PCA start";


    featureVectorz.resize(N_CHANNELS); // 14 matrices

    int n = featureVectorzEnsemble.at(0).rows();// any matrix rows ...
    int nFeatures = featureVectorzEnsemble.at(0).cols();// any matrix cols ...

    qDebug() << "samples " << n << "   features " << nFeatures;

    MatrixXd centered;
    MatrixXd cov;
    SelfAdjointEigenSolver<MatrixXd> eig; // eigen solver declaration
    MatrixXd eigValuesMatrix(nFeatures,N_CHANNELS);

    ensembleMeanz.resize(N_CHANNELS);
    ensembleEigVectorz.resize(N_CHANNELS);

    for(int channel=0 ; channel < N_CHANNELS ; channel++)
    {
        ensembleMeanz[channel]=featureVectorzEnsemble.at(channel).colwise().mean();
        centered = featureVectorzEnsemble.at(channel).rowwise() - ensembleMeanz.at(channel);

        cov = ( centered.adjoint() * centered )/(n-1);

        // remember to try put eigen solver declaration here if you face problems ...
        eig.compute(cov);

        // ...
        qDebug() << eig.eigenvalues().rows();
        qDebug() << eig.eigenvalues().cols();
        // ...
        eigValuesMatrix.col(channel)=eig.eigenvalues().reverse(); // reverse to make order from high to low
        ensembleEigVectorz[channel]=eig.eigenvectors().rowwise().reverse(); // reverse to make order from high to low

    }

    // calculate #PCs (principal components) here
    eigValuezMean=eigValuesMatrix.rowwise().mean();

    qDebug() << "eigValuezMean rows " <<eigValuezMean.rows();
    qDebug() << "eigValuezMean cols " <<eigValuezMean.cols();

    nPCz=1;

    // I use tail as the eigen values are sorted ascending (from low to high)
    while(( eigValuezMean.head(nPCz).cwiseAbs().sum() / eigValuezMean.cwiseAbs().sum() ) < threshold)
    {
        nPCz++;
    }

    // take same #PCs from all matrices

    for(int channel=0 ; channel < N_CHANNELS ; channel++)
    {
        centered = featureVectorzEnsemble.at(channel).rowwise() - ensembleMeanz.at(channel);
        featureVectorz[channel]=centered*ensembleEigVectorz.at(channel).leftCols(nPCz);
    }

    qDebug() << "PCA end";
}

void trainFLDz()
{
      qDebug() << "train FLDz start";

      //qDebug() << "PCs  " << nPCz;

      QVector<MatrixXd> pcFeatureVectorz; // feature vectors based on principal components

      pcFeatureVectorz.resize(nPCz);

      //qDebug() << "---- 1 ----";

      // initialize matices of pcFeatureVectorz

      for(int pc=0 ; pc < nPCz ; pc++)
      {
          pcFeatureVectorz[pc].setZero(nCharz*N_ROWCOL,N_CHANNELS);
      }

      //qDebug() << "---- 2 ----";
      //qDebug() << featureVectorz.at(0).rows();
      //qDebug() << featureVectorz.at(0).cols();


      for(int pc=0 ; pc < nPCz ; pc++)
      {
          for(int channel=0 ; channel < N_CHANNELS; channel++) // pc ... principal components
          {
              pcFeatureVectorz[pc].col(channel)=featureVectorz.at(channel).col(pc);
          }

      }

      //qDebug() << "---- 3 ----";

      // set number of weight vectors
      weightVectorz.resize(nPCz); // remember that weightVectos is QVector of MatrixXd

      for(int pc=0 ; pc < nPCz ; pc++)
      {

       /** FLD Method **/


          MatrixXd XP300; // feature vectors of P300 class
          MatrixXd XNP300; // feature vectors of non-P300 class


          // set sizes
          XP300.resize(nCharz*2,pcFeatureVectorz.at(pc).cols()); // we have 2 out of 12 P300 signals for each character
          XNP300.resize(nCharz*10,pcFeatureVectorz.at(pc).cols()); // we have 10 out of 12 non-P300 signals for each character

          int row_xp300 = 0;
          int row_xnp300 = 0;

          for(int r=0 ; r < pcFeatureVectorz.at(pc).rows() ; r++)
          {
              if(labelz(r)==1)
              {
                  XP300.row(row_xp300)=pcFeatureVectorz.at(pc).row(r);
                  row_xp300++;
              }
              else
              {
                  XNP300.row(row_xnp300)=pcFeatureVectorz.at(pc).row(r);
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


        weightVectorz[pc] = Sw.colPivHouseholderQr().solve((M1-M2).transpose());
      }

      //qDebug() << weightVectorz.at(0).rows() <<"  " << weightVectorz.at(0).cols();

      qDebug() << "train FLDz end";
}


void getTrainingScores()
{
    qDebug() << "Training Scores";


    // initialize training scores
    trainingScores.setZero(nCharz*N_ROWCOL,nPCz); //vectors * scores for each classifier

    // forming feature vectors

    QVector<MatrixXd> pcFeatureVectorz; // feature vectors based on principal components

    pcFeatureVectorz.resize(nPCz);

    //qDebug() << "---- 1 ----";

    // initialize matices of pcFeatureVectorz

    for(int pc=0 ; pc < nPCz ; pc++)
    {
        pcFeatureVectorz[pc].setZero(nCharz*N_ROWCOL,N_CHANNELS);
    }

    //qDebug() << "---- 2 ----";

    for(int pc=0 ; pc < nPCz ; pc++)
    {
        for(int channel=0 ; channel < N_CHANNELS ; channel++) // pc ... peincipal components
        {
            pcFeatureVectorz[pc].col(channel) = featureVectorz.at(channel).col(pc);
        }

    }

    //qDebug() << "---- 3 ----";

    for(int pc=0 ; pc < nPCz ; pc++)
    {
        trainingScores.col(pc) = pcFeatureVectorz[pc]*weightVectorz.at(pc);
    }


    qDebug() << trainingScores.rows() << " " << trainingScores.cols();

    qDebug() << "Training Scores end";
}


void trainFinalFLD()
{
      qDebug() << "train Final FLD start";

      // set size of Weight Vector
      Weight.resize(nPCz);

       /** FLD Method **/

      MatrixXd XP300; // feature vectors of P300 class
      MatrixXd XNP300; // feature vectors of non-P300 class

      // set sizes
      XP300.resize(nCharz*2,trainingScores.cols()); // we have 2 out of 12 P300 signals for each character
      XNP300.resize(nCharz*10,trainingScores.cols()); // we have 10 out of 12 non-P300 signals for each character

      int row_xp300 = 0;
      int row_xnp300 = 0;

      for(int r=0 ; r < trainingScores.rows() ; r++)
      {
          if(labelz(r)==1)
          {
              XP300.row(row_xp300)=trainingScores.row(r);
              row_xp300++;
          }
          else
          {
              XNP300.row(row_xnp300)=trainingScores.row(r);
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


      Weight = Sw.colPivHouseholderQr().solve((M1-M2).transpose());
      //bias=-(M1+M2)*Weight/2; // you can add 1 to features instead of usig bias and to reduce code change
      // remember you drop the bias as it is not needed in max score :)

      //qDebug() << weightVectorz.at(0).rows() <<"  " << weightVectorz.at(0).cols();

      qDebug() << "train Final FLD end";
}




void testFinalFLD()
{
    qDebug() << "test Final FLD start";

    // forming feature vectors

    MatrixXd testingScores; // to be defined


    //qDebug() << "---- 3 ----";

    // testing feature vectors using weights ...
    VectorXd cols;
    VectorXd rows;
    double maxColScore;
    double maxRowScore;
    VectorXd::Index maxRow, maxCol;

    //qDebug() << "---- 4 ----";

    for(int r=0 ; r< nCharz*N_ROWCOL ; r+=N_ROWCOL) // step = 12 (12 vectors / char)
    {
        cols.setZero(N_COLS);
        rows.setZero(N_ROWS);

        //qDebug() << "---- 5 ----";

        cols += testingScores.block(r,0,N_COLS,testingScores.cols())*Weight;//+bias


        rows += testingScores.block(r+N_COLS,0,N_ROWS,testingScores.cols())*Weight;//+bias

        maxColScore =cols.maxCoeff(&maxCol);
        maxRowScore=rows.maxCoeff(&maxRow);

        labelzOfTestedData.append(N_COLS*maxRow+maxCol);

    }

    qDebug() << "test Final FLD end";
}



// can be removed as we just want to train
void testFLDz()
{
    qDebug() << "test FLDz start";

    // forming feature vectors

    QVector<MatrixXd> pcFeatureVectorz; // feature vectors based on principal components

    pcFeatureVectorz.resize(nPCz);

    //qDebug() << "---- 1 ----";

    // initialize matices of pcFeatureVectorz

    for(int pc=0 ; pc < nPCz ; pc++)
    {
        pcFeatureVectorz[pc].setZero(nCharz*N_ROWCOL,N_CHANNELS);
    }

    //qDebug() << "---- 2 ----";

    for(int pc=0 ; pc < nPCz ; pc++)
    {
        for(int channel=0 ; channel < N_CHANNELS ; channel++) // pc ... peincipal components
        {
            pcFeatureVectorz[pc].col(channel) = featureVectorz.at(channel).col(pc);
        }

    }

    //qDebug() << "---- 3 ----";

    // testing feature vectors using weights ...
    VectorXd cols;
    VectorXd rows;
    double maxColScore;
    double maxRowScore;
    VectorXd::Index maxRow, maxCol;

    // take the last nPCz values
    VectorXd PCsWeight=eigValuezMean; // PC weight (i.e. eigen values)

    //qDebug() << "---- 4 ----";

    for(int r=0 ; r< nCharz*N_ROWCOL ; r+=N_ROWCOL) // step = 12 (12 vectors / char)
    {
        //maxColScore=0;
        //maxRowScore=0;
        cols.setZero(N_COLS);
        rows.setZero(N_ROWS);

        //qDebug() << "---- 5 ----";

        for(int pc=0 ; pc < nPCz ; pc++)
        {
            cols += (pcFeatureVectorz[pc].block(r,0,N_COLS,N_CHANNELS)*weightVectorz.at(pc))/PCsWeight.head(pc+1).sum();


            rows += (pcFeatureVectorz[pc].block(r+N_COLS,0,N_ROWS,N_CHANNELS)*weightVectorz.at(pc))/PCsWeight.head(pc+1).sum();

        }

        maxColScore =cols.maxCoeff(&maxCol);
        maxRowScore=rows.maxCoeff(&maxRow);

        labelzOfTestedData.append(N_COLS*maxRow+maxCol);

    }

    qDebug() << "test FLDz end";
}


// can be removed as we just want to train
void printTex()
{
    qDebug() << "print text start";

//    QString str="A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,0,1,2,3,4,5,6,7,8,9";
//    QStringList charzScreen = str.split(',');

//    QTextStream out(stdout);
//    for(int i=0;i<labelzOfTestedData.size();i++)
//    {
//        out << charzScreen.at(labelzOfTestedData.at(i));
//    }

    qDebug() << "print text end";
}


void saveWeightVectorz()
{

    qDebug()<<"ensembleWeightVectors";

    QString path="/sdcard/BrainTexEd/par_linear/ensembleWeightVectors";
    QFile output(path);

    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    // save number of weight vectors (i.e. = nPCz) at first line
    stream << weightVectorz.size() << endl;
    // save size at 2nd line
    stream << weightVectorz.at(0).size() << endl;

    for(int i=0;i<weightVectorz.size();i++)
    {
        for(int j=0; j < weightVectorz.at(i).size() ; j++)
        {
            stream << weightVectorz.at(i)(j) << " ";
        }

        stream << endl;
    }


    output.close();

}


void saveVectorzMeanz()
{

    qDebug()<<"ensembleMeans";

    QString path="/sdcard/BrainTexEd/par_linear/ensembleMeans";
    QFile output(path);
    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    // save size at first line
    stream << ensembleMeanz.at(0).size() << endl; // it is always 128/dec_fator

    for(int channel=0;channel<N_CHANNELS;channel++)
    {
        for(int i=0;i<ensembleMeanz.at(channel).size();i++)
        {
            stream << ensembleMeanz.at(channel)(i) << " ";
        }
        stream << endl;
    }

    output.close();

}

void saveEigValuezMean()
{

    qDebug()<<"ensembleEigValuesMean";

    QString path="/sdcard/BrainTexEd/par_linear/ensembleEigValuesMean";
    QFile output(path);

    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    // save size at first line
    stream << eigValuezMean.size() << endl; // it is always 128/dec_fator

    for(int i=0;i<eigValuezMean.size();i++)
    {
            stream << eigValuezMean(i) << endl;
    }

    output.close();
}

void saveEigVectorz()
{

    qDebug()<<"ensembleEigVectors";

    QString path="/sdcard/BrainTexEd/par_linear/ensembleEigVectors";
    QFile output(path);

    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    // save #PCs at first line
    stream << nPCz << endl;
    // save size at 2nd and 3rd lines
    stream << ensembleEigVectorz.at(0).rows() << endl;
    stream << ensembleEigVectorz.at(0).cols() << endl;

    for(int channel=0;channel<N_CHANNELS;channel++)
    {

        for(int r=0;r<ensembleEigVectorz.at(channel).rows();r++)
        {
            for(int c=0;c<ensembleEigVectorz.at(channel).cols();c++)
            {
                stream << ensembleEigVectorz.at(channel)(r,c) << " ";
            }

            stream << endl;
        }

    }

    output.close();
}




void saveFinalWeight()
{

    qDebug()<< "Weight";

    QString path="/sdcard/BrainTexEd/par_linear/Weight";
    QFile output(path);

    output.open(QFile::WriteOnly);
    QTextStream stream(&output);

    // save size at first line
    stream << Weight.size() << endl;

    for(int i=0;i<Weight.size();i++)
    {
            stream << Weight(i) << endl;
    }

    output.close();
}


void train_linear()
{
    // list all operations in order
    countSamplez();
    getLabelz();
    loadSignalzData();
    applyCARfilterz();
    applyBandPassFilterz();
    loadFlashingVector();
    loadStimulusCodeVector();
    createDataz();
    createEnsembleFeatureVectorz();
    pcaz(THRESHOLD);
    trainFLDz();
    getTrainingScores();
    trainFinalFLD();    
    //testFinalFLD();
    //printTex();
    saveWeightVectorz();
    saveVectorzMeanz();
    saveEigValuezMean();
    saveEigVectorz();
    saveFinalWeight();
}
