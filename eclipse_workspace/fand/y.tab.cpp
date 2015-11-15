#ifndef lint
static char const 
yyrcsid[] = "$FreeBSD: src/usr.bin/yacc/skeleton.c,v 1.28 2000/01/17 02:04:06 bde Exp $";
#endif
#include <stdlib.h>
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define YYLEX yylex()
#define YYEMPTY -1
#define yyclearin (yychar=(YYEMPTY))
#define yyerrok (yyerrflag=0)
#define YYRECOVERING() (yyerrflag!=0)
static int yygrowstack();
#define YYPREFIX "yy"
#line 2 "fand.y"
/* This grammar parses a single Allele enty of the B11 langauge.*/
/* It returns a PMF<Allele> corresponding the the parsed expression*/
/* 'F' is represented by any empty PMF.*/
/* 'D' is represented by Allele::unknown within the PMF, with its corresponding probability*/
/* it is up to the caller to process the 'D'*/

/* TODO: clean up memory in case of error*/

#include <stdio.h>
#include <string.h>

#include "fand/Parser.h"

extern "C"
{
    int yywrap()
    {
        return 1;
    }
}

extern char *yytext;

int yylex(void);  
int yyparse(void);

void yyerror(const char *str)
{
    /*fprintf(stderr,"syntax error: %s at '%s'\n", str, yytext);*/
}

#line 35 "fand.y"
typedef union 
{
    float        decimal;
    PMF<Allele>* pmf;
} YYSTYPE;
#line 55 "y.tab.cpp"
#define YYERRCODE 256
#define ALLELE 257
#define PROB 258
const short yylhs[] = {                                        -1,
    0,    1,    1,    1,    1,    1,    1,    3,    3,    2,
    4,    4,    6,    5,    5,
};
const short yylen[] = {                                         2,
    1,    1,    1,    3,    1,    1,    3,    1,    3,    1,
    1,    3,    5,    1,    3,
};
const short yydefred[] = {                                      0,
   10,    2,    0,    0,    0,    1,    0,    0,    8,    0,
    0,   14,    0,    0,    0,    0,    4,    0,    0,   12,
    9,    7,   15,    0,   13,
};
const short yydgoto[] = {                                       5,
    6,    7,    8,    9,   13,   10,
};
const short yysindex[] = {                                    -40,
    0,    0,  -60, -252,    0,    0,  -58,  -39,    0,  -57,
 -248,    0,  -38, -247, -252, -246,    0, -252,  -51,    0,
    0,    0,    0,  -52,    0,
};
const short yyrindex[] = {                                      0,
    0,    0,   15,    0,    0,    0,    1,   17,    0,   18,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,
};
const short yygindex[] = {                                      0,
    0,   -2,    0,    4,    0,    0,
};
#define YYTABLESIZE 217
const short yytable[] = {                                       4,
   11,   12,   19,   11,    1,   14,   16,   15,   18,   17,
   20,   22,   24,   25,    3,   23,    5,    6,   21,    0,
    0,    0,    0,    0,    0,    0,    0,    3,    0,    2,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,   11,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    1,
};
const short yycheck[] = {                                      40,
    0,    4,   41,   64,  257,   64,   64,   47,   47,  258,
  258,  258,   64,   66,    0,   18,    0,    0,   15,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   68,   -1,   70,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   47,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,  257,
};
#define YYFINAL 5
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 258
#if YYDEBUG
const char * const yyname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,"'('","')'",0,0,0,0,0,"'/'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"'@'",0,
"'B'",0,"'D'",0,"'F'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
"ALLELE","PROB",
};
const char * const yyrule[] = {
"$accept : top",
"top : pmf",
"pmf : 'F'",
"pmf : 'D'",
"pmf : 'D' '@' PROB",
"pmf : wlist",
"pmf : bterm",
"pmf : bterm '@' PROB",
"wlist : wterm",
"wlist : wlist '/' wterm",
"allele : ALLELE",
"wterm : allele",
"wterm : allele '@' PROB",
"bterm : '(' slist ')' '@' 'B'",
"slist : allele",
"slist : slist '/' allele",
};
#endif
#if YYDEBUG
#include <stdio.h>
#endif
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 10000
#define YYMAXDEPTH 10000
#endif
#endif
#define YYINITSTACKSIZE 200
int yydebug;
int yynerrs;
int yyerrflag;
int yychar;
short *yyssp;
YYSTYPE *yyvsp;
YYSTYPE yyval;
YYSTYPE yylval;
short *yyss;
short *yysslim;
YYSTYPE *yyvs;
int yystacksize;
/* allocate initial stack or double stack size, up to YYMAXDEPTH */
static int yygrowstack()
{
    int newsize, i;
    short *newss;
    YYSTYPE *newvs;

    if ((newsize = yystacksize) == 0)
        newsize = YYINITSTACKSIZE;
    else if (newsize >= YYMAXDEPTH)
        return -1;
    else if ((newsize *= 2) > YYMAXDEPTH)
        newsize = YYMAXDEPTH;
    i = yyssp - yyss;
    newss = yyss ? (short *)realloc(yyss, newsize * sizeof *newss) :
      (short *)malloc(newsize * sizeof *newss);
    if (newss == NULL)
        return -1;
    yyss = newss;
    yyssp = newss + i;
    newvs = yyvs ? (YYSTYPE *)realloc(yyvs, newsize * sizeof *newvs) :
      (YYSTYPE *)malloc(newsize * sizeof *newvs);
    if (newvs == NULL)
        return -1;
    yyvs = newvs;
    yyvsp = newvs + i;
    yystacksize = newsize;
    yysslim = yyss + newsize - 1;
    return 0;
}

#define YYABORT goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR goto yyerrlab

#ifndef YYPARSE_PARAM
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG void
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif	/* ANSI-C/C++ */
#else	/* YYPARSE_PARAM */
#ifndef YYPARSE_PARAM_TYPE
#define YYPARSE_PARAM_TYPE void *
#endif
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG YYPARSE_PARAM_TYPE YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL YYPARSE_PARAM_TYPE YYPARSE_PARAM;
#endif	/* ANSI-C/C++ */
#endif	/* ! YYPARSE_PARAM */

int
yyparse (YYPARSE_PARAM_ARG)
    YYPARSE_PARAM_DECL
{
    register int yym, yyn, yystate;
#if YYDEBUG
    register const char *yys;

    if ((yys = getenv("YYDEBUG")))
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = (-1);

    if (yyss == NULL && yygrowstack()) goto yyoverflow;
    yyssp = yyss;
    yyvsp = yyvs;
    *yyssp = yystate = 0;

yyloop:
    if ((yyn = yydefred[yystate])) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yyssp >= yysslim && yygrowstack())
        {
            goto yyoverflow;
        }
        *++yyssp = yystate = yytable[yyn];
        *++yyvsp = yylval;
        yychar = (-1);
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;
#if defined(lint) || defined(__GNUC__)
    goto yynewerror;
#endif
yynewerror:
    yyerror("syntax error");
#if defined(lint) || defined(__GNUC__)
    goto yyerrlab;
#endif
yyerrlab:
    ++yynerrs;
yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yyssp]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yyssp, yytable[yyn]);
#endif
                if (yyssp >= yysslim && yygrowstack())
                {
                    goto yyoverflow;
                }
                *++yyssp = yystate = yytable[yyn];
                *++yyvsp = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yyssp);
#endif
                if (yyssp <= yyss) goto yyabort;
                --yyssp;
                --yyvsp;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = (-1);
        goto yyloop;
    }
yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    yyval = yyvsp[1-yym];
    switch (yyn)
    {
case 1:
#line 53 "fand.y"
{
            /* return results and clean up temporaries*/
            *Parser::m_pmf = *yyvsp[0].pmf;
            Parser::deletePMF( yyvsp[0].pmf );
        }
break;
case 2:
#line 61 "fand.y"
{
            yyval.pmf = Parser::newPMF(); /* empty.*/
        }
break;
case 3:
#line 66 "fand.y"
{
            yyval.pmf = Parser::newPMF();
            yyval.pmf->insert(std::make_pair(Allele(Allele::unknown), 1.0));
        }
break;
case 4:
#line 72 "fand.y"
{
            yyval.pmf = Parser::newPMF();
            yyval.pmf->insert(std::make_pair(Allele(Allele::unknown), yyvsp[0].decimal));
        }
break;
case 5:
#line 78 "fand.y"
{
            yyval.pmf = yyvsp[0].pmf;
        }
break;
case 6:
#line 83 "fand.y"
{
            yyval.pmf = yyvsp[0].pmf;
        }
break;
case 7:
#line 88 "fand.y"
{
            yyval.pmf = yyvsp[-2].pmf;
            *yyval.pmf *= yyvsp[0].decimal;
        }
break;
case 8:
#line 95 "fand.y"
{
            yyval.pmf = yyvsp[0].pmf;
        }
break;
case 9:
#line 100 "fand.y"
{
            yyval.pmf = yyvsp[-2].pmf;
            *yyval.pmf += *yyvsp[0].pmf;
            Parser::deletePMF( yyvsp[0].pmf );
        }
break;
case 10:
#line 108 "fand.y"
{
            int repeats = int(yyvsp[0].decimal);
            int variant = int( (yyvsp[0].decimal - repeats + 0.01) * 10); /* nasty!*/
            Allele a(repeats, variant);

            yyval.pmf = Parser::newPMF();

            /* check if allele is known            */
            double f = Parser::background(a);
            if (f>0)
            {
                yyval.pmf->insert(std::make_pair(a, 1.0));
            }
            /* else unknown allele - ignore*/
        }
break;
case 11:
#line 126 "fand.y"
{
            yyval.pmf = yyvsp[0].pmf;
        }
break;
case 12:
#line 131 "fand.y"
{
            yyval.pmf = yyvsp[-2].pmf;
            *yyval.pmf *= yyvsp[0].decimal;
        }
break;
case 13:
#line 138 "fand.y"
{
            yyval.pmf = yyvsp[-3].pmf;
            PMF<Allele>::iterator it, current;
            for(it = yyval.pmf->begin(); it != yyval.pmf->end();)
            {
                current = it++;
                double f = Parser::background(current->first);
                if (f>0)
                {
                    current->second *= f;
                }
                else
                {
                    /* f = 0: ignore allele*/
                    yyval.pmf->erase(current);
                }
            }
            yyval.pmf->normalize();
        }
break;
case 14:
#line 160 "fand.y"
{
            yyval.pmf = yyvsp[0].pmf;
        }
break;
case 15:
#line 165 "fand.y"
{
            yyval.pmf = yyvsp[-2].pmf;
            *yyval.pmf += *yyvsp[0].pmf;
            Parser::deletePMF( yyvsp[0].pmf );
        }
break;
#line 521 "y.tab.cpp"
    }
    yyssp -= yym;
    yystate = *yyssp;
    yyvsp -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yyssp = YYFINAL;
        *++yyvsp = yyval;
        if (yychar < 0)
        {
            if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yyssp, yystate);
#endif
    if (yyssp >= yysslim && yygrowstack())
    {
        goto yyoverflow;
    }
    *++yyssp = yystate;
    *++yyvsp = yyval;
    goto yyloop;
yyoverflow:
    yyerror("yacc stack overflow");
yyabort:
    return (1);
yyaccept:
    return (0);
}
