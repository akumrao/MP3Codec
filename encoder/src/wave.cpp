/*
 * http://soundfile.sapp.org/doc/WaveFormat/
 * 
 * 
 */

/* 
 * File:   wave.cpp
 * Author: root
 * 
 * Created on March 27, 2018, 7:50 PM
 */


/* Errorcodes
  -------------
    0 = No Error (OK).
   -1 = Unsupported filetype.
   -2	= Couldn't open file.
   -3 = Unexpected end of file.
   -4 = (lame_global_struct file is not in the format the file-extension says).
   -5 = Important chunk missing.
   -6 = Samples are in unsupported (compressed?) format.
 */


#define  STR_COMM  0x4d4d4f43
#define  STR_SSND  0x444e5353
#define  STR_FORM  0x4d524f46
#define  STR_AIFF  0x46464941
#define  STR_RIFF  0x46464952
#define  STR_WAVE  0x45564157
#define  STR_fmt   0x20746d66
#define  STR_data  0x61746164

#include "wave.h"

wave::wave() {
}

wave::wave(const wave& orig) {
}

wave::~wave() {
}



/*____ openInput() ____________________________________________________________*/

int wave::openInput(lame_global_struct * psInfo, char * pFileName){
    int x;
    char header[3 * 4];

   

    /* Set Filepointer */

    if (pFileName == NULL)
    {
        psInfo->fp = stdin;
    } else
    {
        psInfo->fp = fopen(pFileName, "rb");
        if (psInfo->fp == NULL) goto couldNotOpen;
    }


    /* Read and analyze header */

    if (fread(header, 4, 3, psInfo->fp) != 3) goto couldNotOpen;

    if (intlLong(&header[0]) == STR_RIFF && intlLong(&header[8]) == STR_WAVE)
        x = initWAV(psInfo);
    //else if (intlLong(&header[0]) == STR_FORM && intlLong(&header[8]) == STR_AIFF)
       // x = initAIFF(psInfo);
    else
    {
       // memcpy(psInfo->preReadBuffer, header, 12);
       // psInfo->nPreReadBytes = 12;
       // x = initRAW(psInfo);
    }

    if (x == FALSE)
    {
   
        fclose(psInfo->fp);
      
        return FALSE;
    }


    /* Set some flags */

   // if (psInfo->fchannels)
   //     psInfo->outputType = STEREO;
  //  else
  //      psInfo->outputType = DOWNMIX_MONO;

   // psInfo->outputFreq = psInfo->freq;

    return TRUE;

couldNotOpen:
  //  psInfo->errcode = -2;
   // psInfo->samplesLeft = 0;
    return FALSE;
}

/*____ readSamples() __________________________________________________________*/

int wave::readSamples(lame_global_struct * psInfo, uint nSamples, short * wpSamples){
    int retVal = 0;
    uint i;
    uint readSamples;
    char temp;
    short tmp;
/*
    if (psInfo->fchannels == TRUE && (psInfo->outputType == DOWNMIX_MONO
            || psInfo->outputType == LEFT_CHANNEL_MONO || psInfo->outputType == RIGHT_CHANNEL_MONO))
        readSamples = nSamples * 2;
    else
        readSamples = nSamples;


    if (psInfo->samplesLeft == 0)
        return 0;

    if (psInfo->samplesLeft != 0xFFFFFFFF)
    {
        if (readSamples < psInfo->samplesLeft)
            psInfo->samplesLeft -= readSamples;
        else
        {
            readSamples = psInfo->samplesLeft;
            psInfo->samplesLeft = 0;
        }
    }

    if (psInfo->filetype == WAV)
        retVal = readWAVSamples(psInfo, readSamples, wpSamples);
    //else if (psInfo->filetype == AIFF)
        //retVal = readAIFFSamples(psInfo, readSamples, wpSamples);
    else if (psInfo->filetype == RAW)
        retVal = readRAWSamples(psInfo, readSamples, wpSamples);

    if (psInfo->samplesLeft == 0 || retVal == FALSE)
    {
        psInfo->samplesLeft = 0;
        if (psInfo->fp != stdin)
            fclose(psInfo->fp);
    }


    ///Possibly swap byteorder 

    if (psInfo->bits == 16 && psInfo->byteorder != BYTEORDER)
    {
        for (i = 0; i < readSamples; i++)
        {
            temp = ((char *) wpSamples)[i * 2];
            ((char *) wpSamples)[i * 2] = ((char *) wpSamples)[i * 2 + 1];
            ((char *) wpSamples)[i * 2 + 1] = temp;
        }
    }


   //  Convert between 8/16-bit 

    if (psInfo->bits == 8)
    {
        for (i = readSamples - 1; i > 0; i--)
            wpSamples[i] = ((short) ((unsigned char *) wpSamples)[i]) << 8;
        wpSamples[i] = ((short) ((unsigned char *) wpSamples)[i]) << 8;
    }

    //Convert unsigned to signed 

    if (psInfo->fSign == FALSE)
    {
        for (i = 0; i < readSamples; i++)
            wpSamples[i] ^= 0x8000;
    }

    // Convert from Stereo to Mono or inverse stereo in a number of ways 

    if (psInfo->outputType != STEREO && psInfo->fchannels == TRUE)
    {
        if (psInfo->outputType == DOWNMIX_MONO)
            for (i = 0; i < readSamples / 2; i++)
                wpSamples[i] = (short) ((((int) wpSamples[i * 2]) + ((int) wpSamples[i * 2 + 1])) >> 1);

        if (psInfo->outputType == LEFT_CHANNEL_MONO)
            for (i = 0; i < readSamples / 2; i++)
                wpSamples[i] = wpSamples[i * 2];


        if (psInfo->outputType == RIGHT_CHANNEL_MONO)
            for (i = 0; i < readSamples / 2; i++)
                wpSamples[i] = wpSamples[i * 2 + 1];

        if (psInfo->outputType == INVERSE_STEREO)
        {
            for (i = 0; i < readSamples; i += 2)
            {
                tmp = wpSamples[i];
                wpSamples[i] = wpSamples[i + 1];
                wpSamples[i + 1] = tmp;
            }
        }
        retVal /= 2;
    }
 */
    return retVal;
}

/*____ closeInput() ___________________________________________________________*/

int wave:: closeInput(lame_global_struct * psInfo){
   
            fclose(psInfo->fp);
       
        return TRUE;
   
}

/*____ initWAV() ______________________________________________________________*/

 int wave:: initWAV(lame_global_struct * psInfo)
{
    char header[3 * 4];
    int fFmtChunkFound = FALSE;

    struct
    {
        short wFormatTag; /* Format category */
        short wChannels; /* Number of channels */
        int dwSamplesPerSec; /* Sampling rate */
        int dwAvgBytesPerSec; /* For buffer estimation */
        short wBlockAlign; /* Data block size */
        short bitsPerSample; /* Actually a PCM-specific additional byte... */
    } sFmtChunk;

    char aTemp[sizeof ( sFmtChunk)];


    /* Go through the chunks until we have found 'data'. */


    if (fread(header, 4, 2, psInfo->fp) != 2)
        return false;

    while (intlLong(&header[0]) != STR_data)
    {
        if (intlLong(&header[0]) == STR_fmt)
        {
            if (fread(aTemp, sizeof ( sFmtChunk), 1, psInfo->fp) != 1)
                return false;
           fseek(psInfo->fp, intlLong(&header[4]) - sizeof ( sFmtChunk), SEEK_CUR);
            fFmtChunkFound = TRUE;
        } else
            fseek(psInfo->fp, intlLong(&header[4]), SEEK_CUR);

        if (fread(header, 4, 2, psInfo->fp) != 2) 
             return false; 
    }


    /* Fill in sFmtChunk */

    sFmtChunk.wFormatTag = intlShort(aTemp);
    sFmtChunk.wChannels = intlShort(aTemp + 2);
    sFmtChunk.dwSamplesPerSec = intlLong(aTemp + 4);
    sFmtChunk.dwAvgBytesPerSec = intlLong(aTemp + 8);
    sFmtChunk.wBlockAlign = intlShort(aTemp + 12);
    sFmtChunk.bitsPerSample = intlShort(aTemp + 14);


    /* Process the data in sFmtChunk */

    if (fFmtChunkFound != TRUE)
    {
        //psInfo->errcode = -5;
        return FALSE;
    }

    if (sFmtChunk.wFormatTag != 1)
    {
        //psInfo->errcode = -6;
        return FALSE; /* Not a PCM-sample. */
    }

    if (sFmtChunk.wChannels > 2)
    {
        //psInfo->errcode = -6;
        return FALSE; /* More than two channels. */
    }

    psInfo->samplerate_in = sFmtChunk.dwSamplesPerSec;
    psInfo->num_channels = sFmtChunk.wChannels ;
    
    int length = intlLong(&header[4]);
     
   psInfo->num_samples = length / ( psInfo->num_channels  * ((sFmtChunk.bitsPerSample + 7) / 8));
            
            

   
    myPrintf("length=%d.\n",length) ;
   // myPrintf("samplesLeft=%d\n", sFmtChunk.bitsPerSample );
    myPrintf("bits per sample =%d \n", psInfo->num_samples);
    myPrintf("SamplesPerSec =%d\n", psInfo->samplerate_in);
    myPrintf("Channels =%d\n", psInfo->num_channels );
  
    return TRUE;



}

/*____ readWAVsamples() _______________________________________________________*/

 uint wave::readWAVSamples(lame_global_struct * psInfo, int nSamples, short * wpSamples){
    return fread(wpSamples, 16 / 8, nSamples, psInfo->fp);
}

/*____ initRAW() ______________________________________________________________*/

 int wave::initRAW(lame_global_struct * psInfo){

    /* By default we think it is ... */
/*
    psInfo->freq = 44100;
    psInfo->length = 0xFFFFFFFF;
    psInfo->samplesLeft = 0xFFFFFFFF;
    psInfo->fchannels = TRUE;
    psInfo->bits = 16;
    psInfo->fSign = TRUE;
    psInfo->byteorder = BYTEORDER;
    psInfo->filetype = RAW;
*/
    return TRUE;
}

/*____ readRAWsamples() _______________________________________________________*/

 uint wave::readRAWSamples(lame_global_struct * psInfo, int nSamples, short * wpSamples){
    int nPreReadSamples = 0;
/*
    if (psInfo->nPreReadBytes != 0)
    {
        memcpy(wpSamples, psInfo->preReadBuffer, psInfo->nPreReadBytes);
        wpSamples += psInfo->nPreReadBytes / 2;

        nPreReadSamples = psInfo->nPreReadBytes / (psInfo->bits / 8);
        psInfo->nPreReadBytes = 0;
    }
    return fread(wpSamples, psInfo->bits / 8, nSamples - nPreReadSamples, psInfo->fp) + nPreReadSamples;
    */
    
    return 0;
}



/*____ myFseek() ______________________________________________________________*/

/* We can't use the real fseek() since you can't seek in a stream (stdin) */
/*
int myFseek(FILE * fp, int offset){
    char dummy[256];

    while (offset / 256)
    {
        fread(dummy, 256, 1, fp);
        offset -= 256;
    }

    if (offset)
        fread(dummy, offset, 1, fp);

    return 0;
}

*/