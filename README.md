


# Teste Bioinformatica 2024

Por que trabalhar na NeoGenomica ?
===============================

A NeoGenomica entra no mercado como um dos principais laboratórios do Brasil a oferecer tecnologia de sequenciamento genético de nova geração, focada na identificação, análise e diagnóstico de doenças raras. Além disso, disponibilizamos uma ampla gama de exames genéticos voltados para análises clínicas.

Nossa equipe é composta por especialistas renomados em suas áreas, como biomédicos, médicos e bioinformatas, que utilizam as mais avançadas tecnologias e equipamentos para realizar testes genéticos de ponta. Em especial, o time de bioinformática lida diariamente com ferramentas, plataformas e sistemas que manipulam grandes volumes de dados gerados a partir do sequenciamento de DNA e bancos de dados sobre doenças genéticas.

Com o desafio de processar esses grandes volumes de dados, desenvolvemos ferramentas de alto desempenho e especializadas, que apoiam nossos analistas na realização de diagnósticos clínicos precisos. Na NeoGenomica, trabalhamos com Big Data voltado para a saúde, o que reforça nosso compromisso com a inovação. Fazer parte do nosso time significa colaborar com profissionais de diversas formações, unidos pelo objetivo de oferecer a milhões de brasileiros acesso a informações detalhadas e profundas sobre sua saúde, auxiliando-os em suas escolhas de vida.



## Teste técnico para processo seletivo NeoGenomica



Etapa 1
---------
Nesta primeira etapa queremos avaliar seus conhecimentos teóricos sobre o tema de bioinformática e sua capacidade de sintetizar as respostas. Para evitarmos respostas pré-prontas em CHAT GPT, queremos que vocÊ exponha em um material de 4 slides para nós as respostas destas 4 perguntas a seguir. Você deverá nos enviar os slides por e-mail em formato PDF. Deveremos na entrevista de retorno pedir para que você exponha o seu material com as respostas em até 12 min.

1. **O que são variantes genéticas e como elas podem ser classificadas?**
2. **Descreva o pipeline típico para a detecção de variantes em um experimento de sequenciamento de exoma. Quais são as etapas principais?**
3. **Quais ferramentas bioinformáticas você utilizaria para a chamada de variantes em dados de sequenciamento de DNA humano? Explique brevemente a função de cada uma.**
4. **O que é VCF (Variant Call Format) e como ele é utilizado na análise de variantes genéticas?**



Etapa 2
---------

Neste teste você deverá construir um pipeline de bioinformática usando linguagem do seu interesse (shell script, r , python, nextflow) para detecção e anotação de variantes oriundos de dados brutos de NGS DNASeq.  Os dados são de uma amostra de controle humano de sexo feminino. Neste teste você deverá desenvolver o pipeline seguindo as seguintes etapas:

 - Alinhamento das sequências de DNA (FASTQs)
 - Chamada e detecção de variantes SNVs e INDELs
 - Anotação de Variantes
 

Instruções
----------

- Realize o fork deste projeto para que crie um espelho em seu repositório (ex: github.com/marcelcaraciolo/bioinfotest) github. Mais instruções de como fazer o fork [aqui](https://docs.github.com/pt/free-pro-team@latest/github/getting-started-with-github/fork-a-repo).

- Os dados brutos das amostras se encontram on-line será necessário realizar o download das mesmas. Elas estão em formato FASTQ.gz. Instruções de como baixar os fastqs estão na pasta `data` - `data_files_download.sh`

- Coloque todo o código realizado dentro da pasta `code` e os resultados coloque numa pasta `output` (Arquivos BAM, SAM, VCF, Arquivo de respostas).

- Há um questionário de perguntas dentro da pasta `output` com nome `QUESTION.txt` , responda as perguntas dentro do arquivo, salve e commit dentro do seu repositório quando concluído. Estas respostas são obrigatórias e farão parte de sua avaliação técnica.

- Não vamos precisar executar o código aqui localmente do seu pipeline, mas vamos querer ver como você realizou todo o processo desde o alinhamento até a chamada de variantes, portanto fica claro que não serão aceitas soluções em plataformas on-line automatizadas de pipeline como Galaxy, etc. O seu código pode ser colocado em um ou mais arquivos, fica à seu critério de como organizar o código do pipeline.

- Iremos utilizar o genoma de referência da UCSC hg19.fasta para alinhamento e chamada de variantes, para auxiliar o processo deixamos o link aqui dos arquivos necessários para esta etapa.

- Para as etapas de processamento , nossa recomendação de ferramentas são:

  -  BWA (http://bio-bwa.sourceforge.net/) para etapa de alinhamento
  -  FreeBayes (https://github.com/freebayes/freebayes) para etapas de chamada de variantes. Será necessário enviar um parâmetro com o arquivo das regiões-alvo de interesse (``--target``) , para que ele não rode o algoritmo de detecção em todo o genoma humano.  Disponibilizamos o arquivo de regiões neste repositório em ``data``: ``BRCA.list``.
  - snpeff para anotação funcional das variantes (https://pcingola.github.io/SnpEff/)
  
 - Para agilizar o desenvolvimento, criamos um arquivo ``Dockerfile`` a partir do sistema docker (https://www.docker.com/) contendo os arquivos necessários para a prova e já com os programas instalados, se seguir por este caminho poderá economizar um enorme tempo! Para facilitar após instalar o docker criamos o arquivo ``run_test.sh``que ao ser executado no terminal, já realiza o build da imagem do Dockerfile e já faz um ssh com a máquina. Rode a partir do diretório `bioinfotest`.  Para baixar os arquivos de amostra e fasta do genoma humano versao 19, rode dentro da pasta `data`o comando `data_files_download.sh`

  
Resultados Esperados
--------------------

- Vamos precisar que sejam enviados os arquivos: BAM file com os alinhamentos, o BAI file (arquivo de índice), o VCF file (arquivo de variantes) e o arquivo anotado em formato VCF.

- O arquivo QUESTION.txt dentro da pasta ``output`` preenchido com as respostas embaixo de cada quesito. 

- Para facilitar ao terminar o seu teste, commit todo o seu projeto no seu respositório forkeado (bifurcado) e nos envie o link do seu repositório junto a resposta do seu teste admissional.
