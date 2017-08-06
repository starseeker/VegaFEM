/* getopts.c */
/* $Header: H:/ADS/SRC/MSG/RCS/GETOPTS.C 1.3 1995/05/01 14:30:33 Reini Exp $ */

/*********************************************************************/
/*                                                                   */
/*  This Program Written by Paul Edwards.                            */
/*                                                                   */
/*********************************************************************/
/*********************************************************************/
/*                                                                   */
/*  getopts - scan the command line for switches.                    */
/*                                                                   */
/*  This program takes the following parameters:                     */
/*                                                                   */
/*  1) argc (which was given to main)                                */
/*  2) argv (which was given to main)                                */
/*  3) Array of options                                              */
/*                                                                   */
/*  Returns the number of the argument that is next to be processed  */
/*  that wasn't recognised as an option.                             */
/*  Example of use:                                                  */
/*                                                                   */
/*  #include <getopts.h>                                             */
/*  int baud = 2400;                                                 */
/*  char fon[13] = "telix.fon";                                      */
/*  opt_t opttable[] =                                               */
/*  {                                                                */
/*    { "b", OPTINT, &baud },                                        */
/*    { "f", OPTSTR, fon },                                          */
/*    { NULL, 0, NULL }                                              */
/*  }                                                                */
/*  optup = getopts(argc,argv,opttable);                             */
/*                                                                   */
/*  The OPTINT means that an integer is being supplied.  OPTSTR      */
/*  means a string (with no check for overflow).  Also there is      */
/*  OPTBOOL which means it is a switch that is being passed, and an  */
/*  OPTLONG to specify a long.                                       */
/*                                                                   */
/*  This program was inspired by a description of a getargs function */
/*  written by Dr Dobbs Small-C Handbook.  Naturally I didn't get    */
/*  to see the code, otherwise I wouldn't be writing this!           */
/*                                                                   */
/*  This program is dedicated to the public domain.  It would be     */
/*  nice but not necessary if you gave me credit for it.  I would    */
/*  like to thank the members of the International C Conference      */
/*  (in Fidonet) for the help they gave me in writing this.          */
/*                                                                   */
/*  Written 16-Feb-1990.                                             */
/*                                                                   */
/*  Additions by Reini Urban, Graz April 1995:                       */
/*    stringoptions (OPTSTR) as /sstring or as /s string             */
/*      after a string must be whitespace! no "/sstring/a" allowed   */
/*      problems with "/a/b/c"                                       */
/*      a '-' after a bool switch means 0 not next switch            */
/*      valid args: "/a-/b+/sstring /w test"                         */
/*      returns: 3 for "test"                                        */
/*    accepts a '=' after the switch                                 */
/*      eg: /l=string                                                */
/*    new type OPTOPT for options with one or more suboptions        */
/*      eg: for optimization /o with suboptions ailmn...             */
/*      returns the suboptions as string                             */
/*                                                                   */
/*  No warranties!                                                   */
/*
      $Log: GETOPTS.C $
      Revision 1.3  1995/05/01 14:30:33  Reini
      new type OPTSUBOPT added,
      delim '=' between opt and value added, eg: /s=string
*/
/*  Additions by Jernej Barbic, Carnegie Mellon University, 2005:    */
/*    delim ' ' between opt and value added, eg: /s string           */
/*    note: function returns argc upon success                       */

/*********************************************************************/

#include "getopts.h"
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "vegalong.h"

int getopts(int argc, char **argv, opt_t opttable[])
{
  int i = 0;
  //char header[] = "$Id: GETOPTS.C 1.3 1995/05/01 14:30:33 Reini Exp $";
  //char credits[] = "written by Paul Edwards, improved by Reini Urban and Jernej Barbic";

  //skip the first arg
  argv++;
  argc--;
  for (i=1;i<=argc;i++)
  {
    char * arg = *argv;
    //printf("%d %s\n",i, arg);
    if ((*arg != '-') && (*arg != '/')) return (i);
    arg++;

    for (int j=0; opttable[j].sw != NULL;j++) //find options for this arg
    {
      //options are case sensitive!!
      int l = strlen(opttable[j].sw);
      bool bool_ = true;
      
      //printf("%s\n", opttable[j].sw);
      if (strncmp(arg,opttable[j].sw, l) == 0) //found this option
      {
      	char * var = NULL; //pointers to the input data for this option

        if ( ((int)strlen(arg) == l) &&
            ( (i+1 <= argc) && (**(argv+1) != '-') ) )//&&
          //( **(argv+1) != '/') )
        { //copy the next arg: /i int
          argv++; i++;
          var = *argv;
          //printf("next arg\n");
          //printf("**:%s\n",opttable[j].var);
        }
        else
        { //copy the rest: /iint
          if ( ( *(arg+l)=='=' ) || (*(arg+l)==' ') ) // accept '=' or ' ' after the switch
            var = arg + l + 1;
          else
            var = arg + l;
          //printf("$$:%s\n",opttable[j].var);
        }

        // if (*var == '\0')
        //   continue;
        //printf("%d %s\n",i, var);
        switch ((int)opttable[j].opttyp)
        {
        case OPTINT :
          // if ( ( *(arg+l)=='=' ) || (*(arg+l)==' ') ) // accept '=' or ' ' after the switch
          //   *((int *)opttable[j].var) = (int)strtol(arg+l+1,NULL,10);
          // else
          //   *((int *)opttable[j].var) = (int)strtol(arg+l,NULL,10);
          // if (errno == ERANGE)
          //   return (i);
          // break;
          if (*var != '\0')
          {
            *((int *)opttable[j].var) = (int)strtol(var,NULL,10);

            if (errno == ERANGE)
              return (i);
          }
          break;

        case OPTLONG :
          if (*var != '\0')
          {
            *((vegalong *)opttable[j].var) = strtol(var,NULL,10);
            if (errno == ERANGE)
              return (i);
          }
          break;

          // case OPTSUBOPT :
          //   //option with more than 1 char as suboptions:
          //   // /oail or /op
          //   strcpy((char *)opttable[j].var, arg+l);
          //   arg += strlen(arg)-1;
          //   break;
        case OPTSTR :
          if (*var != '\0')
            strcpy((char *)opttable[j].var, var);
          break;

        case OPTBOOL :
          //check a + or - after the sw: "/a- /b+"
          if ( *var == '-' )
            bool_ = false;
          else
            bool_ = true;

          *((bool *)opttable[j].var) = bool_;
          break;
        }
        break;      //break the for
      } // end if (strncmp(arg,opttable[j].sw, l) == 0) 
    } // end for (int j=0; opttable[j].sw != NULL;j++)

    argv++;  //try the next arg, eg: "/a /b"
  }
  return (i);
}

