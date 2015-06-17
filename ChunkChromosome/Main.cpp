////////////////////////////////////////////////////////////////////// 
// ChunkChromosome/Main.cpp 
// (c) 2000-2011 Goncalo Abecasis
// 
// This file is distributed as part of the Goncalo source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile Goncalo.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Friday August 05, 2011
// 
 
#include "StringArray.h"
#include "Parameters.h"
#include "InputFile.h"
#include <string>

using std::string;

string getFileName(const string& s) {

       char sep = '/';

#ifdef _WIN32
          sep = '\\';
#endif

             size_t i = s.rfind(sep, s.length( ));
             if (i != string::npos) {
                return(s.substr(i+1, s.length( ) - i));
             }

             return s;
}


int main(int argc, char ** argv)
   {
   printf("ChunkChromosome - Splits Data File into Chunks to Facilitate Windowed Imputation\n");
   #ifdef VERSION
      printf("This version is stamped " VERSION "\n");
   #endif
   printf("\n");

   String datfile;
   int    chunkSize = 5000;
   int    overlap = 500;

   ParameterList pl;

   pl.Add(new StringParameter('d', "Data File", datfile));
   pl.Add(new IntParameter('n', "Chunk Size", chunkSize));
   pl.Add(new IntParameter('o', "Overlap", overlap));
   pl.Read(argc, argv);
   pl.Status();

   StringArray input;
   input.Read(datfile);

   int markers = 0;

   StringArray tokens;
   for (int i = 0; i < input.Length(); i++)
      {
      tokens.ReplaceTokens(input[i]);

      if (tokens.Length() != 2 || tokens[0].SlowCompare("M") != 0)
         continue;

      markers++;
      }

   printf("Data files includes %d markers ...\n", markers);

   if (markers == 0)
      {
      printf("   *** NOTHING TO DO ***\n\n");
      exit(-1);
      }

   int chunks = markers / chunkSize;

   if (chunks == 0) chunks = 1;

   printf("   Data will be split into %d chunks (~%d markers per chunk)\n\n",
          chunks, (markers + chunks - 1) / chunks);

   String filename;
   String datfilename;
   datfilename.printf("%s",(getFileName(datfile.c_str())).c_str()); 
   filename.printf("autoChunk-%s", (const char *) datfilename);

   IFILE autoChunk = ifopen(filename, "wb");

   ifprintf(autoChunk, "\n"
                       "START\tSTOP\tCORE_START\tCORE_END\n");

   if (autoChunk == NULL)
      {
      printf("   Problem opening file %s\n", (const char *) filename);
      return(-2);
      }

   for (int c = 0; c < chunks; c++)
      {
      int start = int((long long) (c) * markers / chunks);
      int stop  = int((long long) (c + 1) * markers / chunks) - 1;

      int oStart = start - overlap;
      int oStop = stop + overlap;

      if (oStart < 0) oStart = 0;
      if (oStop >= markers) oStop = markers - 1;

      String filename;
      filename.printf("chunk%d-%s", c + 1, (const char *) datfilename);

      IFILE output = ifopen(filename, "wb");
      IFILE snps = ifopen(filename + ".snps", "wb");

      if (output == NULL)
         {
         printf("   Problem opening file %s\n", (const char *) filename);
         return(-2);
         }

      if (snps == NULL)
         {
         printf("   Problem opening file %s.snps\n", (const char *) filename);
         return(-2);
         }

      int marker = 0;

      String mStart("start"), mStop("stop");
      String mFirst, mLast;

      for (int i = 0; i < input.Length(); i++)
         {
         tokens.ReplaceTokens(input[i]);

         if (tokens.Length() != 2 || tokens[0].SlowCompare("M") != 0)
            {
            ifprintf(output, "%s\n", (const char *) input[i]);
            continue;
            }

         if (marker < oStart || marker > oStop)
            ifprintf(output, "S2 %s\n", (const char *) tokens[1]);
         else
            {
            ifprintf(output, "%s\n", (const char *) input[i]);
            ifprintf(snps, "%s\n", (const char *) tokens[1]);
            }

         if (marker == oStart)
            mStart = tokens[1];

         if (marker == start)
            mFirst = tokens[1];

         if (marker == stop)
            mLast = tokens[1];

         if (marker == oStop && oStop < markers)
            mStop = tokens[1];

         marker++;
         }

      ifclose(output);

      if (mStart == "start") mStart = mFirst;
      if (mStop == "stop") mStop = mLast;

      ifprintf(autoChunk, "%s\t%s\t%s\t%s\n",
                  (const char *) mStart, (const char *) mStop,
                  start == 0 ? "start" : (const char *) mFirst,
                  stop == markers -1 ? "stop" : (const char *) mLast);

      printf("   Chunk %d written to %s ...\n"
             "       Flanking markers from %s to %s ...\n"
             "       Core region from %s to %s ...\n",
             c + 1, (const char *) filename,
             (const char *) mStart, (const char *) mStop,
             (const char *) mFirst, (const char *) mLast);
      }

   ifclose(autoChunk);

   printf("\nChuncking information stored in file %s ...\n", (const char *) filename);
   } 
