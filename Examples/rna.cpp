/*
 **********************************************************************************
 *  Copyright (C) 2010  Massachusetts Institute of Technology
 *  Copyright (C) 2010  Yuan Tang <yuantang@csail.mit.edu>
 * 		                Charles E. Leiserson <cel@mit.edu>
 * 	 
 *  Written by: Rezaul Alam Chowdhury <rezaul@mit.edu>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

/* Test bench - RNA secondary structure prediction with simple pseudoknots */
/* ( quadratic space implementation of Akutsu's algorithm from Chowdhury et al., TCBB, 2010 ) */

#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
#define SIMPLE 0
#define N_RANK 2

#ifndef base_pair
  #undef base_pair
#endif  

#define base_pair( a1, a2 ) ( a1 + a2 == 0 )

#define match_pair( a1, a2 ) ( base_pair( a1, a2 ) ? 1 : -INF )

#ifndef MAX_SEQ_LENGTH
  #define MAX_SEQ_LENGTH  1000000
#endif

#ifndef INF
  #define INF ( 10 * MAX_SEQ_LENGTH )
#endif

#define pArrayR2T3  Pochoir_Array< int, N_RANK, 3 >
#define P_ARRAY_R2_T3  pArrayR2T3

//#define TEST_TYPEDEF

//#ifdef TEST_TYPEDEF
//  #define P_ARRAY_R2_T3 pArrayR2T3
//#else
//  #define P_ARRAY_R2_T3 Pochoir_Array< int, N_RANK, 3 >
//#endif  

enum alphabet{ A = -1, U = 1, G = -2, C = 2, ALPHABET_SIZE };

enum err_msgs{ SEQUENCE_READ, LENGTH_READ, NO_SEQUENCE, SEQUENCE_TOO_LONG, INVALID_SEQUENCE_LENGTH, FILE_OPEN_ERROR, MEM_ALLOC_ERROR };


int read_next_seq( char *fn, char **sq, int *len )
{
  FILE *fp;

  if ( fn != NULL )
    {
     if ( ( fp = fopen( fn, "rt" ) ) == NULL ) return FILE_OPEN_ERROR;
    }
  else fp = stdin;

  if ( fscanf( fp, "%d", len ) != 1 ) return NO_SEQUENCE;
  
  ( *sq ) = ( char * ) malloc( ( *len + 2 ) * sizeof( char ) );
  
  if ( *sq == NULL ) return MEM_ALLOC_ERROR;

  int i = 1, c;
  
  while ( ( i <= *len ) && ( ( c = fgetc( fp ) ) != NULL ) )
    {
      ( *sq )[ i++ ] = c;
    } 
  
  ( *sq )[ i ] = 0;

  if ( fn != NULL ) fclose( fp );

  return SEQUENCE_READ;
}


int read_sequence( char *fnX, int *nx, char **X )
{
  int l;

  l = read_next_seq( fnX, X, nx );

  if ( l == FILE_OPEN_ERROR )
    {
      printf( "Error: Unable to open input file ( %s )!\n\n", fnX );
      return 0;
    }
  
  if ( l == MEM_ALLOC_ERROR )
    {
      printf( "Error: Memory allocation error!\n\n" );
      return 0;
    }
  
  if ( l == NO_SEQUENCE )
    {
      printf( "Error: Failed to read sequence!\n\n" );
      return 0;
    }

  return 1;
}



int initRandSeq( int nX, char **X )
{
    ( *X ) = ( char * ) malloc( ( nX + 2 ) * sizeof( char ) );
  
    if ( ( *X == NULL ) ) 
      {
        printf( "Error: Memory allocation error!\n\n" );      
        return 0;
      }   

    char *SYM = "AUGC";    
    int nS = 4;//sizeof( SYM ) / sizeof( SYM[ 0 ] );

    srand( time( NULL ) );

    for ( int i = 1; i <= nX; ++i )
       ( *X )[ i ] = SYM[ rand( ) % nS ];

    ( *X )[ nX + 1 ] = 0;       
    
    return 1;
}


void convertToInts( int nX, char *X )
{
  for ( int i = 1; i <= nX; ++i )
    {
      switch ( X[ i ] )
        {
          case 'A' : X[ i ] = A; break;

          case 'U' : X[ i ] = U; break;

          case 'G' : X[ i ] = G; break;

          case 'C' : X[ i ] = C; break;
        }

    }
}



void print_usage( char *prog )
{
  printf( "Usage: %s [ options ]\n\n", prog );

  printf( "Options:\n" );

  printf( "\t-s fname : input file containing the seqeunece (file format: sequence length in line 1 followed by the sequence in line 2)\n" );
  printf( "\t-r len   : generate a random sequence of length len\n\n" );

  printf( "\t-d       : run standard DP algorithm\n" );
  printf( "\t-i       : run iterative stencil\n\n" );  
   
  printf( "\t-h       : print this help screen\n\n" );
}



int read_command_line( int argc, char *argv[ ], int &nx, char **fnX, int &run_standard_dp, int &run_iter_stencil )
{
  nx = 0;
  run_standard_dp = 0;
  run_iter_stencil = 0;

  for ( int i = 1; i < argc; )
    {
     int j = i;

     if ( !strcasecmp( argv[ i ], "-d" ) )
       {
        run_standard_dp = 1;
         
        i++;
        
        if ( i >= argc ) break;
       }
     

     if ( !strcasecmp( argv[ i ], "-i" ) )
       {
        run_iter_stencil = 1;
         
        i++;
        
        if ( i >= argc ) break;
       }

     
     if ( !strcasecmp( argv[ i ], "-r" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           printf( "Error: Missing sequence length ( specify -r sequence-length )!\n\n" );
           return 0;
          }

        nx = atoi( argv[ i + 1 ] );

        if ( nx < 1 )
          {
           printf( "Error: Specify positive sequence length!\n\n" );
           return 0;
          }
          
        if ( nx > MAX_SEQ_LENGTH )
          {
           printf( "Error: Cannot generate sequences longer than %d!\n\n", MAX_SEQ_LENGTH );
           return 0;
          }          

        i += 2;

        if ( i >= argc ) break;
       }
       
     if ( !strcasecmp( argv[ i ], "-s" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           printf( "Error: Missing input file name ( specify -s sequence-filename )!\n\n" );
           return 0;
          }

        ( *fnX ) = strdup( argv[ i + 1 ] );

        i += 2;
        
        if ( i >= argc ) break;
       }
                     
     if ( !strcasecmp( argv[ i ], "-h" ) || !strcasecmp( argv[ i ], "-help" ) || !strcasecmp( argv[ i ], "--help" ) )
       {
        print_usage( argv[ 0 ] );
        exit( 0 );
        i++;

        if ( i >= argc ) break;
       }


     if ( i == j )
       {
        printf( "Error: Unknown option ( %s )!\n\n", argv[ i ] );
        return 0;
       }
    }
   
   return 1;
}


void stencilRNAi0( int nX, char *X, int i_0, 
                   P_ARRAY_R2_T3 &SL, 
                   P_ARRAY_R2_T3 &SR, 
                   P_ARRAY_R2_T3 &SM,
                   P_ARRAY_R2_T3 &SMAX,                    
                   P_ARRAY_R2_T3 &SP )
{
    Pochoir< int, N_RANK, 3 > pRNA;    
//    Pochoir_Domain I( 0, nX + 1 ), K( 0, nX + 1 );
    Pochoir_Shape< N_RANK > pRNA_shape[ 6 ] = { { 2, 0, 0 }, 
                                                { 1, 0, 0 }, { 1, 0, -1 }, { 1, -1, 0 },
                                                { 0, -1, 0 }, { 0, 0, -1 } };    

    for ( int k_0 = 1; k_0 <= nX; ++k_0 )
      for ( int i = i_0; i < k_0 - 1; ++i )
         SP( 0, i, k_0 ) = -INF;

    for ( int t = 0; t < 2; ++t )
      for ( int i = 0; i <= nX; ++i )
        for ( int k = 0; k <= nX; ++k )      
          {
            SMAX( t, i, k ) = -INF;
            
            int j = t - i - k, jj = nX - j + 1;
            
            if ( ( j >= 0 ) && ( j <= nX ) && ( i_0 - 1 <= i ) && ( i < jj ) && ( jj <= k ) )
               {                   
                 if ( jj == k )
                   {
                     if ( base_pair( X[ i ], X[ jj ] ) ) SL( t, i, k ) = 1;
                     else SL( t, i, k ) = -INF;             
                   }
                 else if ( i >= i_0 )
                        {
                          if ( base_pair( X[ i ], X[ jj ] ) ) SL( t, i, k ) = 1;
                          else SL( t, i, k ) = -INF;
                        }  
                      else SL( t, i, k ) = 0;
                         
                 if ( i == i_0 - 1 )
                   {
                     if ( k == jj + 1 )
                       {
                         if ( base_pair( X[ jj ], X[ k ] ) ) SR( t, i, k ) = 1;
                         else SR( t, i, k ) = -INF;                                             
                       }
                     else SR( t, i, k ) = 0;  
                   }
                 else   
                   {
                     if ( base_pair( X[ jj ], X[ k ] ) ) SR( t, i, k ) = 1;
                     else SR( t, i, k ) = -INF;             
                   }
        
                 if ( i == i_0 - 1 ) SM( t, i, k ) = 0;  
                 else if ( t > 0 )
                         {
                           int v1 = -INF, v2 = -INF;
                           
                           if ( jj < k ) 
                             {
                               v1 = max( SR( 0, i, k - 1 ), SM( 0, i, k - 1 ) );
                               v1 = max( v1, SMAX( 0, i, k ) );
                             } 
                           
                           if ( i - 1 < jj ) v2 = max( SL( 0, i - 1, k ), SM( 0, i - 1, k ) );                           
                                                      
                           SM( t, i, k ) = max( v1, v2 );  
                         }
                      else SM( t, i, k ) = 0;     
                                     
                 int v = max( SL( t, i, k ), SR( t, i, k ) );       
                        
                 SMAX( t, i, k ) = max( v, SM( t, i, k ) );       
                            
                 if ( ( i_0 <= i ) && ( jj < k ) )           
                   {
                     SP( 0, i, k ) = max( SMAX( t, i, k ), SP( 0, i, k ) );
                   }
                 else SP( 0, i, k ) = -INF;   
               }  
          }


    Pochoir_kernel_2D( pRNA_fn, t, i, k )
      
       int j = t + 2 - i - k, jj = nX - j + 1;

       if ( ( j >= 0 ) && ( j <= nX ) && ( i_0 - 1 <= i ) && ( i < jj ) && ( jj <= k ) )
         {                   
           if ( jj == k )
             {
               if ( base_pair( X[ i ], X[ jj ] ) ) SL( t + 2, i, k ) = 1;
               else SL( t + 2, i, k ) = -INF;             
             }
           else if ( i >= i_0 )
                  {
                    if ( base_pair( X[ i ], X[ jj ] ) ) 
                      {
                        if ( i - 1 < jj + 1 )  SL( t + 2, i, k ) = 1 + SMAX( t, i - 1, k );
                        else SL( t + 2, i, k ) = 1;
                      }  
                    else SL( t + 2, i, k ) = -INF;
                 }  
               else SL( t + 2, i, k ) = 0;
                   
           if ( i == i_0 - 1 )
             {
               if ( k == jj + 1 )
                 {
                   if ( base_pair( X[ jj ], X[ k ] ) ) SR( t + 2, i, k ) = 1;
                   else SR( t + 2, i, k ) = -INF;                                             
                 }
               else SR( t + 2, i, k ) = 0;  
             }
           else   
             {
               if ( base_pair( X[ jj ], X[ k ] ) ) 
                 {
                   if ( jj + 1 <= k - 1 ) SR( t + 2, i, k ) = 1 + SMAX( t, i, k - 1 );
                   else SR( t + 2, i, k ) = 1;
                 }  
               else SR( t + 2, i, k ) = -INF;             
             }

           if ( i == i_0 - 1 ) SM( t + 2, i, k ) = 0;  
           else
             {
               int v1 = -INF, v2 = -INF;
               
               if ( jj < k ) 
                 {
                   v1 = max( SR( t + 1, i, k - 1 ), SM( t + 1, i, k - 1 ) );
                   v1 = max( v1, SMAX( t + 1, i, k ) );
                 } 
               
               if ( i - 1 < jj ) v2 = max( SL( t + 1, i - 1, k ), SM( t + 1, i - 1, k ) );
               
               SM( t + 2, i, k ) = max( v1, v2 );  
             }
                               
           int v = max( SL( t + 2, i, k ), SR( t + 2, i, k ) );       
                  
           SMAX( t + 2, i, k ) = max( v, SM( t + 2, i, k ) );       
                      
           if ( ( i_0 <= i ) && ( jj < k ) ) SP( 0, i, k ) = max( SMAX( t + 2, i, k ), SP( 0, i, k ) );
           else SP( 0, i, k ) = -INF;   
         }
                      	      
    Pochoir_kernel_end

    pRNA.registerShape( pRNA_shape );

    pRNA.registerArray( SL );
    pRNA.registerArray( SR );
    pRNA.registerArray( SM );
    pRNA.registerArray( SMAX );
    pRNA.registerArray( SP );
                    
//    pRNA.registerDomain( I, K );    

    int t = 3 * nX - 1;

    pRNA.run( t, pRNA_fn );
}




void iterativeStencilRNAi0( int nX, char *X, int i_0, 
                   P_ARRAY_R2_T3 &SL, 
                   P_ARRAY_R2_T3 &SR, 
                   P_ARRAY_R2_T3 &SM, 
                   P_ARRAY_R2_T3 &SMAX,                    
                   P_ARRAY_R2_T3 &SP )
{
    for ( int k_0 = 1; k_0 <= nX; ++k_0 )
      for ( int i = i_0; i < k_0 - 1; ++i )
         SP( 0, i, k_0 ) = -INF;

    for ( int t = 0; t < 2; ++t )
      for ( int i = 0; i <= nX; ++i )
        for ( int k = 0; k <= nX; ++k )      
          {
            SMAX( t, i, k ) = -INF;
                      
            int j = t - i - k, jj = nX - j + 1;
            
            if ( ( j >= 0 ) && ( j <= nX ) && ( i_0 - 1 <= i ) && ( i < jj ) && ( jj <= k ) )
               {                   
                 if ( jj == k )
                   {
                     if ( base_pair( X[ i ], X[ jj ] ) ) SL( t, i, k ) = 1;
                     else SL( t, i, k ) = -INF;             
                   }
                 else if ( i >= i_0 )
                        {
                          if ( base_pair( X[ i ], X[ jj ] ) ) 
                            {
                              if ( i - 1 < jj + 1 )  SL( t + 2, i, k ) = 1 + SMAX( t, i - 1, k );
                              else SL( t + 2, i, k ) = 1;
                            }  
                          else SL( t + 2, i, k ) = -INF;
                        }  
                      else SL( t, i, k ) = 0;
                         
                 if ( i == i_0 - 1 )
                   {
                     if ( k == jj + 1 )
                       {
                         if ( base_pair( X[ jj ], X[ k ] ) ) SR( t, i, k ) = 1;
                         else SR( t, i, k ) = -INF;                                             
                       }
                     else SR( t, i, k ) = 0;  
                   }
                 else   
                   {
                     if ( base_pair( X[ jj ], X[ k ] ) ) SR( t, i, k ) = 1;
                     else SR( t, i, k ) = -INF;             
                   }
        
                 if ( i == i_0 - 1 ) SM( t, i, k ) = 0;  
                 else if ( t > 0 )
                         {
                           int v1 = -INF, v2 = -INF;

                           if ( jj < k ) 
                             {
                               v1 = max( SR( 0, i, k - 1 ), SM( 0, i, k - 1 ) );
                               v1 = max( v1, SMAX( 0, i, k ) );
                             } 
                           
                           if ( i - 1 < jj ) v2 = max( SL( 0, i - 1, k ), SM( 0, i - 1, k ) );
                                          
                           SM( t, i, k ) = max( v1, v2 );  
                         }
                      else SM( t, i, k ) = 0;     
                                     
                 int v = max( SL( t, i, k ), SR( t, i, k ) );       
                        
                 SMAX( t, i, k ) = max( v, SM( t, i, k ) );       
                            
                 if ( ( i_0 <= i ) && ( jj < k ) )           
                   {
                     SP( 0, i, k ) = max( SMAX( t, i, k ), SP( 0, i, k ) );
                   }
                 else SP( 0, i, k ) = -INF;   
               }  
          }

    
    for ( int t = 0; t <= 3 * nX - 2; ++t )
      cilk_for ( int i = 0; i < nX + 1; ++i )
        for ( int k = 0; k < nX + 1; ++k )      
          {
             int j = t + 2 - i - k, jj = nX - j + 1;
      
             if ( ( j >= 0 ) && ( j <= nX ) && ( i_0 - 1 <= i ) && ( i < jj ) && ( jj <= k ) )
               {                   
                 if ( jj == k )
                   {
                     if ( base_pair( X[ i ], X[ jj ] ) ) SL( t + 2, i, k ) = 1;
                     else SL( t + 2, i, k ) = -INF;             
                   }
                 else if ( i >= i_0 )
                        {
                          if ( base_pair( X[ i ], X[ jj ] ) ) SL( t + 2, i, k ) = 1 + SMAX( t, i - 1, k );
                          else SL( t + 2, i, k ) = -INF;
                       }  
                     else SL( t + 2, i, k ) = 0;
                         
                 if ( i == i_0 - 1 )
                   {
                     if ( k == jj + 1 )
                       {
                         if ( base_pair( X[ jj ], X[ k ] ) ) SR( t + 2, i, k ) = 1;
                         else SR( t + 2, i, k ) = -INF;                                             
                       }
                     else SR( t + 2, i, k ) = 0;  
                   }
                 else   
                   {
                     if ( base_pair( X[ jj ], X[ k ] ) ) 
                       {
                         if ( jj + 1 <= k - 1 ) SR( t + 2, i, k ) = 1 + SMAX( t, i, k - 1 );
                         else SR( t + 2, i, k ) = 1;
                       }  
                     else SR( t + 2, i, k ) = -INF;             
                   }
        
                 if ( i == i_0 - 1 ) SM( t + 2, i, k ) = 0;  
                 else
                   {
                     int v1 = -INF, v2 = -INF;

                     if ( jj < k ) 
                       {
                         v1 = max( SR( t + 1, i, k - 1 ), SM( t + 1, i, k - 1 ) );
                         v1 = max( v1, SMAX( t + 1, i, k ) );
                       } 
                     
                     if ( i - 1 < jj ) v2 = max( SL( t + 1, i - 1, k ), SM( t + 1, i - 1, k ) );
                                                         
                     SM( t + 2, i, k ) = max( v1, v2 );  
                   }
                                     
                 int v = max( SL( t + 2, i, k ), SR( t + 2, i, k ) );       
                        
                 SMAX( t + 2, i, k ) = max( v, SM( t + 2, i, k ) );       
                            
                 if ( ( i_0 <= i ) && ( jj < k ) ) SP( 0, i, k ) = max( SMAX( t + 2, i, k ), SP( 0, i, k ) );
                 else SP( 0, i, k ) = -INF;   
               }
          }
}


int stencilRNA( int nX, char *X, bool recursive )
{
    P_ARRAY_R2_T3 SL( nX + 1, nX + 1 ), SR( nX + 1, nX + 1 ), SM( nX + 1, nX + 1 ), SMAX( nX + 1, nX + 1 );
    P_ARRAY_R2_T3 SP( nX + 1, nX + 1 ), S( nX + 2, nX + 2 );

    for ( int i_0 = 1; i_0 <= nX; ++i_0 )
      {
        if ( recursive ) stencilRNAi0( nX, X, i_0, SL, SR, SM, SMAX, SP );
        else iterativeStencilRNAi0( nX, X, i_0, SL, SR, SM, SMAX, SP );

        for ( int k_0 = 1, v = -INF; k_0 <= nX; ++k_0 )
          {
            v = max( SP( 0, i_0, k_0 ), v );
        
            for ( int i = i_0 + 1; i < k_0 - 1; ++i )
               v = max( SP( 0, i, k_0 ), v );
        
            SP( 0, i_0, k_0 ) = v; 
          } 
        
        for ( int k_0 = 0; k_0 < i_0 + 2; ++k_0 )
           SP( 1, i_0, k_0 ) = -INF;
                
        for ( int k_0 = i_0 + 2; k_0 <= nX; ++k_0 )
           SP( 1, i_0, k_0 ) = max( SP( 1, i_0, k_0 - 1 ), SP( 0, i_0, k_0 ) );            
      }
      
    for ( int j = 0; j <= nX + 1; ++j )
       S( 0, nX + 1, j ) = 0;

    for ( int i = 0; i <= nX + 1; ++i )
       S( 0, i, 0 ) = 0;

    for ( int i = nX; i >= 1; --i )
      for ( int j = i + 1; j <= nX; ++j )
        {
          int v;
          
          if ( base_pair( X[ i ], X[ j ] ) ) v = 1 + S( 0, i + 1, j - 1 );
          else v = -INF;
          
          v = max( v, SP( 1, i, j ) );
          
          for ( int k = i + 1; k <= j; k++ )
            v = max( v, S( 0, i, k - 1 ) + S( 0, k, j ) );
            
          S( 0, i, j ) = v; 
        }  

    return S.interior( 0, 1, nX );
}



int main( int argc, char *argv[ ] )
{
    printf( "\nStencil-based DP for RNA secondary structure prediction with simple pseudoknots ( run with option -h for help ).\n\n" );

    int nX;
    char *X = NULL, *RNA = NULL;    
    char *fnX = NULL;
    int runStandardDP, runIterativeStencil;

    if ( !read_command_line( argc, argv, nX, &fnX, runStandardDP, runIterativeStencil ) )
      {
        print_usage( argv[ 0 ] );
        return 1;
      }

    int gotSeq = 0;

    if ( nX < 1 )
      {
       if ( fnX == NULL )
         {
           printf( "Error: Missing input file name ( specify -s sequence-filename )!\n\n" );
           print_usage( argv[ 0 ] );
           return 1;         
         }
          
       gotSeq = read_sequence( fnX, &nX, &RNA );
      }
    else gotSeq = initRandSeq( nX, &RNA );
        
    RNA[ 0 ] = ' ';
    X = strdup( ( const char * ) RNA );
    convertToInts( nX, X );

    if ( gotSeq )
      {
        printf( "Sequence length = %d\n\n", nX );
      
        struct timeval start, end;

        printf( "Running pochoir-based DP..." );
        fflush( stdout );
               
        gettimeofday( &start, 0 );        
        int maxNumBP = stencilRNA( nX, X, true );    
        gettimeofday( &end, 0 );

        double t0 = tdiff( &end, &start );
              
        printf( "\n\nPochoir:\n" );
        if ( maxNumBP == -INF ) printf( "\t maximum number of base pairs = -inf\n" );    
        else printf( "\t maximum number of base pairs = %d\n", maxNumBP );    
        printf( "\t running time = %.3lf sec\n\n", t0 );    
      
        if ( runIterativeStencil )
          {
            printf( "Running iterative stencil..." );
            fflush( stdout );
                          
            gettimeofday( &start, 0 );        
            int maxNumBPITST = stencilRNA( nX, X, false );
            gettimeofday( &end, 0 );

            double t1 = tdiff( &end, &start );
          
            printf( "\n\nIterative Stencil:\n" );
            if ( maxNumBPITST == -INF ) printf( "\t maximum number of base pairs = -inf\n" );    
            else printf( "\t maximum number of base pairs = %d\n", maxNumBPITST );                
            if ( t0 > 0 ) printf( "\t running time = %.3lf sec ( %.3lf x Pochoir )\n\n", t1, t1 / t0 );    
            else printf( "\t running time = %.3lf sec\n\n", t1 );    
          }
      }

    if ( RNA != NULL ) free( RNA );
    if ( X != NULL ) free( X );
       
    return 0;
}