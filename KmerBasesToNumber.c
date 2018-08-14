#include <string.h>
#include <stdlib.h>
#include <stdio.h>

unsigned long long basesToNumber(char *s)
{
    unsigned long long ret = 0;
    int i = 0;
    while (s[i] != '\0')
    {
        ret = ret << 2;
        switch(s[i])
        {
            case 'A':
                ret += 0;
                break;
            case 'C':
                ret += 1;
                break;
            case 'G':
                ret += 2;
                break;
            case 'T':
                ret += 3;
                break;
        }
        i++;

    }
    return ret;
}

int main(int argc, const char *argv[])
{
    char line[1000];
    char kmer[50], count[100];
    int i;
    while (fgets(line, 999, stdin))
    {
        for (i=0; line[i] != ' '; i++);
        strncpy(kmer, line, i);
        printf("%llu%s",  basesToNumber(kmer), &line[i]);
    }

}
