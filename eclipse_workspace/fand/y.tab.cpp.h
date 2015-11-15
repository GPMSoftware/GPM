#ifndef YYERRCODE
#define YYERRCODE 256
#endif

#define ALLELE 257
#define PROB 258
typedef union 
{
    float        decimal;
    PMF<Allele>* pmf;
} YYSTYPE;
extern YYSTYPE yylval;
