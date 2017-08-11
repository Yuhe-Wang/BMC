#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <string.h>
#include <random>


int main(int argc, char *argv[])
{
	char str[200];
	char chs[100];
	printf("This program will help you create an encrypted file about your gmail account.\n");
	printf("Unfortunately only gmail is supported now.\n\n");
	printf("Input your gmail address: ");
	scanf("%s", str);
	size_t len1 = strlen(str);
	printf("\nInput the password: ");
	scanf("%s", chs);
	printf("\n");
	size_t len2 = strlen(chs);
	strcpy(str + len1 + 1, chs);
	if (len1 > 100 || len2 > 100)
	{
		printf("The length of gmail account name or password should not be over 100 characters!\n");
		getchar();
		return 0;
	}
	
	//fill the rest char with random char
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, 127);
	for (size_t n = len1 + len2 + 2; n < 200; ++n) str[n] = dis(gen);

	//encrypt the str
	const char* key = "simpleXORcipher";
	int keylen = sizeof(key) / sizeof(char);
	for (int i = 0; i < 200; ++i)
	{
		for (int j = 0; j < keylen; ++j)
		{
			str[i] = str[i] ^ key[j];
		}
	}

	FILE* fp = fopen("myEmail.encrypt", "wb");
	if (NULL == fp)
	{
		printf("Error: Cannot write the encrypted file!\n");
		getchar();
		return -1;
	}
	fwrite(str, sizeof(char)*200, 1, fp);
	fclose(fp);
	printf("Success: the encrypted file named myEmail.encrypt was created!\n\n");
	printf("Press any key to continue...\n");
	getchar();
	getchar();
	return 0;
}