%{
#include <stdio.h>
#include <string.h>

#include "fand/PMF.h"
#include "fand/Allele.h"

#include "y.tab.cpp.h"

static YY_BUFFER_STATE bytes_buf;

void start_scan_bytes(char const *buf, size_t size)
{
    bytes_buf = yy_scan_bytes(buf, size);
}

void end_scan_bytes()
{
    yy_delete_buffer(bytes_buf);
}


%}

%%
X                       yylval.decimal=Allele::X; return ALLELE;
Y                       yylval.decimal=Allele::Y; return ALLELE;
[1-9][0-9]*(\.[0-4])?   yylval.decimal=atof(yytext); return ALLELE;
0?\.[0-9]+              yylval.decimal=atof(yytext); return PROB;
[ \t]+                  /* ignore whitespace */;
\n                      return 0; /* logical EOF */
.                       return yytext[0];
%%
