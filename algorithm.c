#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define SEQ_RECORD_MEMORY_CHUNK 40
#define SEQ_SEQUENCE_MEMORY_CHUNK 409600
#define SEQ_MICROSATELLITE_MEMORY_CHUNK 6400

int *readConfig(FILE *f){
    char x[8];
    int *output = malloc(6 * sizeof(int));
    int cursor=0;
    while (fscanf(f, "%1023s", x) == 1){
        if(cursor%2==1)
            output[cursor/2]=atoi(x);
        cursor++;
    }
    int a = 0;
    return output;
}


// Sequence structure

typedef struct {
    char *array;
    size_t used;
    size_t size;
} sequence;

void initSequence(sequence *a, size_t initialSize){
    a->array = malloc(initialSize * sizeof(char));
    a->used = 0;
    a->size = 0;
}

void insertSequence(sequence *a, char element){
    if (a->used == a->size){
        a->size += SEQ_SEQUENCE_MEMORY_CHUNK;
        a->array = realloc(a->array, a->size * sizeof(char));
    }
    a->array[a->used++] = element;
}

// Record structure

typedef struct {
    char *name;
    char *description;
    sequence *sequence;
} record;

typedef struct {
    record *array;
    size_t used;
    size_t size;
} recordArray;

void initRecordArray(recordArray *a, size_t initialSize){
    a = malloc(sizeof(recordArray));
    a->array = malloc(initialSize * sizeof(record));
    a->used = 0;
    a->size = 0;
}

void insertRecordArray(recordArray *a, record element){
    if (a->used == a->size){
        a->size += SEQ_RECORD_MEMORY_CHUNK;
        a->array = realloc(a->array, a->size * sizeof(record));
    }
    a->array[a->used++] = element;
}

// Microsatellite structure

typedef struct {
    char *sequence;
    char *motif;
    int period;
    int repeat;
    int start;
    int end;
    int length;
} microsatellite;

typedef struct {
    microsatellite *array;
    size_t used;
    size_t size;
} microsatelliteArray;

void initMicrosatelliteArray(microsatelliteArray *a, size_t initialSize){
    a->array = malloc(initialSize * sizeof(microsatellite));
    a->used = 0;
    a->size = initialSize;
}

void insertMicrosatelliteArray(microsatelliteArray *a, microsatellite *element){
    if (a->used == a->size){
        a->size += SEQ_MICROSATELLITE_MEMORY_CHUNK;
        a->array = realloc(a->array, a->size*sizeof(microsatellite));
    }
    a->array[a->used++] = *element;
}

void search_perfect_microsatellites(microsatelliteArray *output, record *record, int *minRepeats){
    char *seq = record->sequence->array;
    size_t len;
    int start;
    int length;
    int repeat;
    int i;
    int j;
    char motif[7] = "\0";

    len = record->sequence->used;
    for (i=0; i<len; i++){
        if (seq[i] == 78)
            continue;
        if (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C' && seq[i] != 'G') {
            char *n = malloc(1000000);
            strncpy(n, seq+i, 100);
            continue;
        }

        for (j=1; j<=6; j++) {
            start = i;
            length = j;
            while (start + length < len && seq[i] == seq[i + j] && seq[i] != 78) {
                i++;
                length++;
            }
            repeat = length / j;
            if (repeat >= minRepeats[j - 1]) {
                microsatellite *m = malloc(sizeof(microsatellite));
                m->motif=malloc(16);
                strncpy(m->motif, seq+start, j);
                m->sequence=malloc(64*sizeof(char));
                strcpy(m->sequence, record->name);

                m->period=j;
                m->repeat=repeat;
                m->start=start+1;
                m->end=start+length;
                m->length=length;
                insertMicrosatelliteArray(output, m);
                free(m);

                i=start+length;
                j=0;
            }else{
                i=start;
            }
        }
    }
    int a = 0;
}

void readFastaFile(recordArray *output, FILE *f){
    initRecordArray(output, 5);

    record *tmpRecord;

    char *line = NULL;
    size_t len = 0;
    size_t read;
    while (read = getline(&line, &len, f) != -1){
        if(strlen(line) == 0)
            continue;
        strtok(line, "\n");
        if(line[0] == '>' && output->used != 0){
            printf("ERROR, unexpected case, please check here with that genome");
        }
        if(line[0] == '>'){
            // Sequence name and description
            tmpRecord = malloc(sizeof(record));
            char *x = malloc(64);
            sscanf(line, "%63s", x);
            tmpRecord->name = malloc(64*sizeof(char));
            tmpRecord->description = malloc(1024*sizeof(char));
            strcpy(tmpRecord->name, x);
            strcpy(tmpRecord->description, line+strlen(x));
            tmpRecord->sequence = malloc(sizeof(sequence));
            initSequence(tmpRecord->sequence, SEQ_SEQUENCE_MEMORY_CHUNK);
            int a = 0;
        }else{
            // Fasta reads
            for(int i = 0; line[i] != '\0'; i++)
                insertSequence(tmpRecord->sequence, line[i]);
            printf("");
        }
        printf("");
    }
    // Append tmp value to output
    insertRecordArray(output, *tmpRecord);
}

int main(int argc, char **argv){
    // Parse input
    const char *infile = NULL;
    const char *outfile = NULL;
    const char *configfile = NULL;

    for(int i = 1; i < argc; i++){
        if(strcmp(argv[i-1], "-i") == 0)
            infile = argv[i];
        if(strcmp(argv[i-1], "-o") == 0)
            outfile = argv[i];
        if(strcmp(argv[i-1], "-c") == 0)
            configfile = argv[i];
    }

    // Open the configuration file
    int *minRepeats;
    if(configfile != NULL){
        FILE *fptr;
        fptr = fopen(configfile, "r");
        minRepeats = readConfig(fptr);
    }

    microsatelliteArray *microsatellites = malloc(sizeof(microsatelliteArray));
    initMicrosatelliteArray(microsatellites, SEQ_MICROSATELLITE_MEMORY_CHUNK);
    if(infile != NULL){
        FILE *fptr;
        fptr = fopen(infile, "r");

        // Read fasta file (get records)
        recordArray *records = malloc(sizeof(recordArray));
        readFastaFile(records, fptr);
        search_perfect_microsatellites(microsatellites, records->array, minRepeats);

        // We are here
        int a = 0;
    }
    return 0;
}