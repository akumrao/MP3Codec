/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   wave.h
 * Author: root
 *
 * Created on March 27, 2018, 7:50 PM
 */

#ifndef WAVE_H
#define WAVE_H

#include <stdio.h>
#include <string.h>
#include "system.h"
#include  "dump.h"
#include "mp3.h"

//typedef		unsigned int	uint;

class wave {
public:
    wave();
    wave(const wave& orig);
    virtual ~wave();

    
public:
    
    typedef enum FileType
{
	WAV,
	AIFF,
	RAW
} FileType; 

typedef enum	SampleType
{
	STEREO,
	INVERSE_STEREO,
	DOWNMIX_MONO,
	LEFT_CHANNEL_MONO,
	RIGHT_CHANNEL_MONO
} SampleType;




public:
int	openInput( lame_global_struct * psInfo, char	* pFileName );
int	readSamples(lame_global_struct * psInfo, unsigned int nSamples, short * wpSamples );
int	closeInput( lame_global_struct * psInfo );

private:
    
    /*____ intlLong() _____________________________________________________________*/

inline uint intlLong(char iLong[4]){
    return ((uint) ((uchar*) iLong)[0]) + (((uint) ((uchar*) iLong)[1]) << 8)
            + (((uint) ((uchar*) iLong)[2]) << 16) + (((uint) ((uchar*) iLong)[3]) << 24);
}

/*____ mcLong() _______________________________________________________________*/

inline uint mcLong(char mcLong[4]){
    return ((uint) ((uchar*) mcLong)[3]) + (((uint) ((uchar*) mcLong)[2]) << 8)
            + (((uint) ((uchar*) mcLong)[1]) << 16) + (((uint) ((uchar*) mcLong)[0]) << 24);
}

/*____ intlShort() ____________________________________________________________*/


inline ushort intlShort(char iShort[2]){
    return ((ushort) ((uchar*) iShort)[0]) + (((ushort) ((uchar*) iShort)[1]) << 8);
}

/*____ mcShort() ______________________________________________________________*/

inline ushort mcShort(char mcShort[2]){
    return ((ushort) ((uchar*) mcShort)[1]) + (((ushort) ((uchar*) mcShort)[0]) << 8);
}

private:
     int initWAV(lame_global_struct * psInfo);
 uint readWAVSamples(lame_global_struct * psInfo, int nSamples, short * wpSamples);

 int initRAW(lame_global_struct * psInfo);
 uint readRAWSamples(lame_global_struct * psInfo, int nSamples, short * wpSamples);



 int myFseek(FILE * fp, int offset);




};

#endif /* WAVE_H */

