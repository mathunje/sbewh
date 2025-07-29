#include "parse.h"

void strtrim(char * str)
{
    char * p = str, * start = str;
    while(isspace((unsigned char)*p)) p++;
    if ( p != str ){
        while( *p )
            *str++ = *p++;
        * str = 0;
    } else {
        while ( *p )
            p++;
        str = p;
    }
    if( *start == 0 )
        return;
    p = str - 1;
    while(p > start && isspace((unsigned char)*p)) p--;
    p[1] = 0;
}
