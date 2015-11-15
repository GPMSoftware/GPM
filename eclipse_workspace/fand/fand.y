%{
// This grammar parses a single Allele enty of the B11 langauge.
// It returns a PMF<Allele> corresponding the the parsed expression
// 'F' is represented by any empty PMF.
// 'D' is represented by Allele::unknown within the PMF, with its corresponding probability
// it is up to the caller to process the 'D'

// TODO: clean up memory in case of error

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
    //fprintf(stderr,"syntax error: %s at '%s'\n", str, yytext);
}

%}

%union 
{
    float        decimal;
    PMF<Allele>* pmf;
}

%token <decimal> ALLELE
%token <decimal> PROB
%type <pmf> pmf;
%type <pmf> allele;
%type <pmf> wlist;
%type <pmf> wterm;
%type <pmf> slist;
%type <pmf> bterm;

%%

top   : pmf
        {
            // return results and clean up temporaries
            *Parser::m_pmf = *$1;
            Parser::deletePMF( $1 );
        }
      ;
        
pmf   : 'F'
        {
            $$ = Parser::newPMF(); // empty.
        }

      | 'D'
        {
            $$ = Parser::newPMF();
            $$->insert(std::make_pair(Allele(Allele::unknown), 1.0));
        }

      | 'D' '@' PROB
        {
            $$ = Parser::newPMF();
            $$->insert(std::make_pair(Allele(Allele::unknown), $3));
        }

      | wlist
        {
            $$ = $1;
        }

      | bterm
        {
            $$ = $1;
        }

      | bterm '@' PROB
        {
            $$ = $1;
            *$$ *= $3;
        }
      ;

wlist : wterm
        {
            $$ = $1;
        }
        
      | wlist '/' wterm
        {
            $$ = $1;
            *$$ += *$3;
            Parser::deletePMF( $3 );
        }
      ;
        
allele : ALLELE
        {
            int repeats = int($1);
            int variant = int( ($1 - repeats + 0.01) * 10); // nasty!
            Allele a(repeats, variant);

            $$ = Parser::newPMF();

            // check if allele is known            
            double f = Parser::background(a);
            if (f>0)
            {
                $$->insert(std::make_pair(a, 1.0));
            }
            // else unknown allele - ignore
        }
       ;
        
wterm : allele
        {
            $$ = $1;
        }
        
      | allele '@' PROB
        {
            $$ = $1;
            *$$ *= $3;
        }
      ;

bterm : '(' slist ')' '@' 'B'
        {
            $$ = $2;
            PMF<Allele>::iterator it, current;
            for(it = $$->begin(); it != $$->end();)
            {
                current = it++;
                double f = Parser::background(current->first);
                if (f>0)
                {
                    current->second *= f;
                }
                else
                {
                    // f = 0: ignore allele
                    $$->erase(current);
                }
            }
            $$->normalize();
        }
      ;

slist : allele
        {
            $$ = $1;
        }
        
      | slist '/' allele
        {
            $$ = $1;
            *$$ += *$3;
            Parser::deletePMF( $3 );
        }
      ;
        
      