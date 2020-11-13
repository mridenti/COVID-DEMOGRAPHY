#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <math.h>
#include <stdbool.h>

#ifndef NEA
#define NEA 16
#endif

#define NUM_FAIXAS NEA

int main(int argc, char** argv) {
    FILE *csv;
    FILE *output;
    char filename[512];
	char *input_path = "input\\";
	errno_t err;

    if (argc == 2)
    {
        sprintf_s(filename, "%scenarios\\%s\\initial.csv", input_path, argv[1]);
		if (fopen_s(&csv, filename, "r") == 1)
		{
			perror("Abertura de arquivo de entrada initial.csv falhou. \n");
			return 1;
		}
    }
    else if (argc == 1)
    {
		sprintf_s(filename, "%scenarios\\default\\initial.csv", input_path);
		if (fopen_s(&csv, filename, "r") != 0)
		{
			perror("Abertura de arquivo de entrada initial.csv falhou. \n");
			return 1;
		}
    }
    else {
        printf_s("Uso: \n");
        printf_s("\tcsv_to_input [nome do cenário] \n");
        printf_s("\n");
        printf_s("Exemplos: \n");
        printf_s("1) Usando o cenário \"Default\":\n");
        printf_s("\tcsv_to_input \n");
        printf_s("2) Usando um cenário chamado de \"HardLockdown\":\n");
        printf_s("\tcsv_to_input HardLockdown\n");
        printf_s("\n");
        printf_s("Cenário: \n");
        printf_s("Um cenário é uma pasta dentro de input/cenarios com o nome do cenário e contendo os arquivos .csv referentes aos parâmetros da simulação que serão utilizados");
        printf_s("\t\n");
        return 0;
    }

    sprintf_s(filename, "%s\\generated-input.txt", input_path);
    if (fopen_s(&output,filename, "w") != 0)
	{
		perror("Abertura de arquivo de saida generated-input.txt falhou. \n");
		return 1;
	}

    fprintf_s(output, "# ---------------------------------------------------------- \n");
    fprintf_s(output, "#            INITIAL VALUES PER AGE GROUP                    \n");
    fprintf_s(output, "# ---------------------------------------------------------- \n");
    fprintf_s(output, "\n");

    char line[4096];
    while (fgets(line, 4096, csv) != NULL)
    {
        char* tmp = _strdup(line);

        const char* tok;
		char* posn;
        for (tok = strtok_s(line, ",", &posn); tok && *tok; tok = strtok_s(NULL, ",\n", &posn))
        {
            fprintf_s(output, "%-15s\t", tok);
        }
        fprintf_s(output, "\n");
        free(tmp);
    }
    fclose(csv);

	// Open PARAMETERS file
    if (argc == 2)
    {
        sprintf_s(filename, "%scenarios\\%s\\parameters.csv", input_path, argv[1]);
        printf_s("Lendo: %s\n", filename);
        if(fopen_s(&csv, filename, "r") != 0)
		{
			perror("Abertura de arquivo de parametros paramenters.csv falhou. \n");
			return 1;
		}
    }
    else if(argc == 1)
    {
        sprintf_s(filename, "%scenarios\\default\\parameters.csv", input_path);
		printf_s("Lendo: %s\n", filename);
		if((err = fopen_s(&csv, filename, "r")) != 0)
		{
			perror("Abertura de arquivo de parametros parameters.csv falhou. \n");
			return 1;
		}
    }
    else {
        printf_s("Uso: \n");
        printf_s("csv_to_input [nome do cenário] \n");
        printf_s("\t ou \n");
        printf_s("csv_to_input \n");
        printf_s("\t (Neste caso usará o cenário default) \n");
    }

	fprintf_s(output, "\n\n");
	fprintf_s(output, "# ---------------------------------------------------------- \n");
	fprintf_s(output, "#               PARAMETERS (fixed in time)                   \n");
	fprintf_s(output, "# ---------------------------------------------------------- \n");
	fprintf_s(output, "\n");

	fgets(line, 4096, csv); // Ignorar header
	while (fgets(line, 4096, csv) != NULL)
	{
		char* tmp = _strdup(line);

		const char* tok;
		char* posn;
		for (tok = strtok_s(line, ",", &posn); tok && *tok; tok = strtok_s(NULL, ",\n", &posn))
		{
			fprintf_s(output, "%-15s\t", tok);
		}
		fprintf_s(output, "\n");
		free(tmp);
	}
	fprintf_s(output, "\n");
	fclose(csv);

	// Print gama, gama_QI, gama_QA and beta by day
	FILE *csv_beta_gama;

    if (argc == 2)
    {
        sprintf_s(filename, "%scenarios\\%s\\beta_gama.csv", input_path,  argv[1]);
        printf_s("Lendo: %s\n", filename);
		if (fopen_s(&csv_beta_gama, filename, "r") != 0)
		{
			perror("Abertura de arquivo de parametros beta_gama.csv falhou. \n");
			return 1;
		}
    }
    else if (argc == 1)
    {
        sprintf_s(filename, "%scenarios\\default\\beta_gama.csv", input_path);
		printf_s("Lendo: %s\n", filename);
		if (fopen_s(&csv_beta_gama, filename, "r") != 0)
		{
			perror("Abertura de arquivo de parametros beta_gama.csv falhou. \n");
			return 1;
		}
		sprintf_s(filename, "%scenarios\\default\\beta_gama.csv", input_path);
    }
    else {
        printf_s("Uso: \n");
        printf_s("csv_to_input [nome do cenário] \n");
        printf_s("\t ou \n");
        printf_s("csv_to_input \n");
        printf_s("\t (Neste caso usará o cenário default) \n");
    }

	// Print gama, gama_QI, gama_QA and beta by day
    fgets(line, 4096, csv_beta_gama); // Ignorar header
	int day = -1;

    while (fgets(line, 4096, csv_beta_gama) != NULL)
    {
		char* beta_tmp = _strdup(line);
		bool first_tok = true;
        const char* tok;
		char *posn;
		for (tok = strtok_s(line, ",", &posn); tok && *tok; tok = strtok_s(NULL, ",\n", &posn))
		{
			if (first_tok) {
				if (atoi(tok) > day) {
					day = atoi(tok);
					fprintf_s(output, "# ---------------------------------------------------------- \n");
					fprintf_s(output, "DAY %s \n", tok);

					// Print GAMA
					fprintf_s(output, "GAMA \t\t");
					tok = strtok_s(NULL, ",\n", &posn); // Step
					for (int ii = 0; ii < NEA; ii++) {
						fprintf_s(output, "%-18s ", tok);
						tok = strtok_s(NULL, ",\n", &posn); // Step
					}
					fprintf_s(output, "\n");

					// Print xI
					fprintf_s(output, "xI \t\t\t");
					for (int ii = 0; ii < NEA; ii++) {
						fprintf_s(output, "%-18s ", tok);
						tok = strtok_s(NULL, ",\n", &posn); // Step
					}
					fprintf_s(output, "\n");

					// Print xA
					fprintf_s(output, "xA \t\t\t");
					for (int ii = 0; ii < NEA; ii++) {
						fprintf_s(output, "%-18s ", tok);
						if (ii != (NEA - 1)) tok = strtok_s(NULL, ",\n", &posn); // Step
					}
					fprintf_s(output, "\n");

					// Print BETA
					fprintf_s(output, "BETA \t\t");
				}
				else if (tok[0] != 10) {
					// Padding
					fprintf(output, "     	    ");
				}
				first_tok = false;
			}
			else {
				fprintf(output, "%-18s ", tok);
			}	
		}
		first_tok = true;
		fprintf_s(output, "\n");
		free(beta_tmp);
    }
    fclose(csv_beta_gama);


    return 0;
}
