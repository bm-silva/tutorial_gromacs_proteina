# Tutorial de Dinâmica Molecular Clássica utilizando o pacote GROMACS

## Introdução

Neste tutorial, iremos utilizar um sistema contendo uma lisozima de *Gallus gallus* em solução aquosa e concentração de 0,150 mol/L de NaCl. Esta proteína está identificada no Protein Data Bank (PDB) pelo código 1AKI. 
Este tutorial se baseou em 2 práticas, já aplicadas:

* http://www.mdtutorials.com/gmx/lysozyme/
* Aula prática do Minicurso de Dinâmica Molecular Básica I Curso de Verão em Bioinformática ministrado pela Dra. Deborah Antunes.

Como o modelo a ser utilzizado, é uma estrutura do PDB, é possível que ele possa conter átomos de solventes, ligantes ou outras moléculas. No caso do PDB 1AKI, se encontra moléculas de água. Tenha em mente que é recomendado retirar todos os átomos relacionados ao solvente. Você pode utilziar programas como o Chimera e o PyMOL para realizar esta tarefa!

## Preparação

Vamos converter o arquivo **.pdb** em arquivos do GROMACS (**.gro**). Para isso, vamos utilizar o programa **pdb2gmx**. Este programa tem diversas opções que podem ser vistas através do comando `pdb2gmx –h`. Digite no terminal:

```
gmx pdb2gmx -f model.pdb -o model_processado.gro -ff amber99sb-ildn -water tip3p
```

* **-f**: arquivo de entrada contendo as coordenadas no formato PDB;
* **-p**: arquivo de saída contendo a topologia no formato TOP;
* **-ff**: indicamos qual campo de força vamos utilizar;
* **-water**: indicamos qual o modelo de água será utilizado;
* **-o**: arquivo de saída contendo as modificações feitas.

Foi gerado, também, o **posre.itp**, arquivo que irá aplicar, na topologia da molécula, forças de restrições às moléculas de água, a serem adicionadas daqui a pouco, e da proteína.
Um arquivo bastante importante é o **topol.top**. Nele é possível encontrar algumas informações importantes:

```
; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"

[ moleculetype ]
; Name            nrexcl
Protein_chain_A     3
```

Estas linhas chamam os parâmetros do campo de força **amber99sb-ildn**, indicando que todos os parâmetros ao longo do arquivo são derivados dele, e define o nome da molécula em `name`. Mais informações sobre exclusões podem ser encontradas no manual GROMACS.

A seguir vemos a região `atoms`:

```
[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
; residue   1 LYS rtp NLYS q +2.0
     1         N3      1    LYS      N      1     0.0966      14.01
     2          H      1    LYS     H1      2     0.2165      1.008
     3          H      1    LYS     H2      3     0.2165      1.008
     4          H      1    LYS     H3      4     0.2165      1.008
```

* **nr**: número do átomo;
* **tipo**: tipo do átomo;
* **resnr**: número do resíduo de aminoácido;
* **residue**:nome do resíduo de aminoácido;
* **atom**: nome do átomo;
* **cgnr**: número do grupo de cargas. Os grupos de carga definem unidades de carga inteira que ajudam a acelerar os cálculos
* **charge**: carga do átomo. O **qtot** é uma quantificação total da carga na molécula
* **mass**: massa do átomo;
* **typeB**, **chargeB**, **massB**: Usado para perturbação de energia livre.

Em seguida é encontrado as regiões de `bonds`, `pairs`, `angles` e `dihedrals`. Essas regiões utilizam o número do átomo para poder informar a geometria espacial da molécula: quem está ligado com quem, qual o ângulo diedro formado e assim por diante.

Uma região importante do arquivo é a região das restrições de posição.

```
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif
```

O arquivo **posre.itp** gerado  define uma constante de força usada para manter os átomos no lugar durante as fazes de equílibrio.

Por fim, vemos campos relacionados a topologia da água a ser utilizada e adicionada, logo em seguida a descrição de susas forças de restrição, a topologia dos íons, nome do systema e a descrição do que está contido no sistema atual.

```
; Include water topology
#include "oplsaa.ff/spce.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include generic topology for ions
#include "oplsaa.ff/ions.itp"

[ system ]
; Name
LYSOZYME

[ molecules ]
; Compound        #mols
Protein_A           1
```

## Caixa de Simulação e Solvatação

Nesta etapa iremos criar o ambiente que irá acomodar a nossa molécula e as moléculas de água a serem adicionadas. Para isso, utilizaremos o programa **editconf**:

```
gmx editconf -f model_processado.gro -o model_caixa.gro -c -d 1 -bt cubic
```

* **-c**: centraliza a molécula na caixa;
* **-d**: distância mínima da molécula para os lados da caixa;
* **-bt**: formato da caixa de simulação;
* **-o**: arquivo de saída contendo as modificações feitas;

Com a caixa de simulação criada, vamos agora adicionar o solvente, que neste caso será água, utilizando o programa **solvate**. O comando **solvate** cria uma caixa d’água baseado nas especificações que se passou pelo **editconf**. O  **solvate** irá adicionar o número correto de moléculas de água necessárias para solvatar a caixa nas dimensões especificadas no passo anterior:

```
gmx solvate -cp model_caixa.gro -cs spc216.gro -o model_solvatado.gro -p topol.top
```

* **-cp**: arquivo de configuração da proteína, resultado do passo anterior;
* **-cd**: configuração do solvente, que neste caso é o spc216.gro, parte do pacote GROMACS, que é um modelo genérico utilizado para o modelo de água que vamos utilizar.
* **-o**: arquivo de saída contendo as modificações feitas;
* **-p**: indicamos qual o arquivo de topologia a ser escrito.

Após a solvatação, no final do arquivo **topol.top** foi adicionado agora uma nova molécula referente ao solvente:

```
[ molecules ]
; Compound  #mols
Protein_A       1 
SOL         10832 
```

Vamos agora neutralizar o sistema e adicionar uma concentração salina ao solvente. Para isto, utilizaremos o programa **genion**. Isso é feito através da substituição de átomos do sistema (tanto do soluto como do solvente) por íons (Cloro ou Sódio, por exemplo). Neste caso, substituiremos moléculas do solvente (SOL) pelos íons de acordo com a carga. 
Antes de prosseguirmos, precisamos utilizar o programa **grompp** (GROMACS PreProcessor). Este programa requer um arquivo **.mdp**, onde estão os parâmetros da simulação a ser realizada pelo **mdrun**, ou pelo **genion**: 

```
gmx grompp -f ions.mdp -c model_solvatado.gro -p topol.top -o ions.tpr -maxwarn 1
```

* **-f**: indicamos o arquivo contendo os parâmetros a serem computados e aplicados pelo programa **genion** no próximo comando;
* **-c**: arquivo com a configuração atual do sistema, originado do passo anterior;
* **-p**: indicamos o arquivo de topologia a ser utilizado e modificado;
* **-o**: arquivo de saída que será utilizado no próximo passo contendo as configurações aplicadas nesta etapa.
* **-maxwarn**: indicamos que o aviso gerado não irá atrapalhar esta etapa. Esse aviso é relacionado a sugestão de neutralização do sistema, que será feito no próximo comando.

Agora vamo para o **genion** que irá de fato adicionar tudo que foi pedido ao sistema:

```
echo SOL | gmx genion -s ions.tpr -p topol.top -o model_salgado.gro -pname NA -nname CL -neutral -conc 0.150
```

Para evitarmos interagir com o programa, nós temos que indicar em qual parte da topologia nós iremos substituir os átomos para a adição dos íons. Neste caso, o comando **echo SOL** passa para o comando depois do *|* o que iríamos selecionar.
* **-s**: o arquivo originado do grompp, chamado também de “arquivo de estado (state file)”;
* **-o**: arquivo de saída contendo as modificações feitas;
* **-p**: arquivo de topologia a ser escrito e modificado;
* **-pname**: indicamos qual o íon positivo iremos utilizar;
* **-nname**: indicamos o íon negativo;
* **-neutral**: esse argumento permite a neutralização de cargas do sistemas;
* **-conc**: a concentração molar salina que desejamos.

Perceba, novamente, que no arquivo **topol.top** foram adicionadas novas moléculas:

```
[ molecules ]
; Compound        #mols
Protein_chain_A    1
SOL                10574
NA                 31
CL                 39
```

**Atenção**: o arquivo **.tpr** é pré-compilado com informações de coordenadas, topologia e parâmetros da simulação. Como é um arquivo binário, não pode ser lido num editor de texto, para poder ter acesso ao seu conteúdo é necessário digitar o comando: `gmxdump -s ions.tpr`.

## Minimização de Energia

Vamos agora realizar a minimização de energia do sistema para remoção de maus contatos de van der Waals oriundos das moléculas de água introduzidas e pela rede de ligações de hidrogênio rompidas. Isso pode gerar energias extremas, as quais são liberadas no início da simulação de dinâmica molecular, podendo comprometer a integridade da estrutura inicial. Este passo também requer o uso prévio do “grompp”, pois assim estaremos estabelecendo as condições nas quais a minimização de energia vai ocorrer e também o tipo de minimização utilizado. Neste tutorial, iremos fazer a minimização do tipo Steepest descent. O arquivo de parâmetros encontra-se no diretório (minim.mdp). Digite no terminal:

```
gmx grompp -f minim.mdp -c model_salgado.gro -p topol.top -o em.tpr
```

* **-f**: indicamos o arquivo contendo os parâmetros a serem computados e aplicados pelo programa “genion” no próximo comando;
* **-c**: arquivo com a configuração atual do sistema, originado do passo anterior;
* **-p**: indicamos o arquivo de topologia a ser utilizado e modificado;
* **-o**: arquivo de saída que será utilizado no próximo passo contendo as configurações aplicadas nesta etapa.

Logo em seguida, vamos rodar a minimização utilizando o mdrun:

```
gmx mdrun -v -s em.tpr -deffnm em
```

* **-v**: “verbose mode”, escreve na tela os passos da minimização;
* **-s**: arquivo pré compilado com informações de coordenadas, topologia e parâmetros da simulação. Esse arquivo é adquirido na etapa anterior com o **grompp**;
* **-deffnm**: coloca os nomes dos arquivos de saída com o nomes indicado.

Ao final da minimização de energia teremos os seguintes arquivos: 
* **em.log**: arquivo log (arquivo texto-ASCII do processo de minimização de energia);
* **em.edr**: arquivo (binário) de energia; 
* **em.trr**: arquivo (binário) da trajetória;
* **em.gro**: arquivo final de coordenadas após a EM. 

O parágrafo a seguir é adaptado de: http://www.mdtutorials.com/gmx/lysozyme/05_EM.html
> Existem dois fatores muito importantes a serem avaliados para determinar se a minimização de energia foi bem-sucedida: A primeira é a **energia potencial** (impressa no final do processo de minimização). A **E<sub>pot</sub>** deve ser negativa e (para uma proteína simples em água) na ordem de 10<sup>5</sup>-10<sup>6</sup>, dependendo do tamanho do sistema e do número de moléculas de água. A segunda característica importante é a força máxima, **F<sub>max</sub>**, cujo alvo foi definido em **minim.mdp** (`emtol = 1000,0`) indicando uma **F<sub>max</sub>** alvo não superior a 1000 kJ mol<sup>-1</sup> nm<sup>-1</sup>. É possível chegar a um Epot razoável com **F<sub>max</sub>** > **emtol**. Se isso acontecer, seu sistema pode não ser estável o suficiente para simulação. Avalie por que isso pode estar acontecendo e talvez altere seus parâmetros de minimização (integrador, emstep, etc).

A minimização de energia assegurou que temos uma estrutura de partida razoável, em termos de geometria e orientação do solvente. Para começar uma dinâmica de produção confiável, é importante equilibrar o solvente e os íons em torno da proteína.
O solvente precisa estar orientado corretamente sobre o soluto na temperatura que se pretende estabelecer a simulação. Depois que conseguirmos chegar na temperatura correta, vamos aplicar pressão ao sistema até atingir a densidade adequada.

## Equilibração do Sistema

A etapa de equilíbrio geralmente é realizada em duas etapas. A primeira fase é realizada sob o ensaio NVT (Número constante de partículas, Volume e Temperatura).. O arquivo “nvt.mdp” se encontra na sua pasta de trabalho. Nós usaremos novamente o “grompp” e o “mdrun” assim como fizemos na etapa de minimização de energia:

```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```

* **-f**: indicamos o arquivo contendo os parâmetros a serem computados e aplicados pelo programa **genion** no próximo comando;
* **-c** e **-r**: arquivo com a configuração atual do sistema, originado do passo anterior;
* **-p**: indicamos o arquivo de topologia a ser utilizado e modificado;
* **-o**: arquivo de saída que será utilizado no próximo passo contendo as configurações aplicadas nesta etapa.

```
gmx mdrun -v -deffnm nvt
```

* **-v**: *verbose mode*, escreve na tela os passos da minimização;
* **-s**: arquivo pré compilado com informações de coordenadas, topologia e parâmetros da simulação. Esse arquivo é adquirido na etapa anterior com o *grompp*;
* **-deffnm**: coloca os nomes dos arquivos de saída com o nomes indicado.

Repare que os arquivos de saída são similares ao da etapa de minimização, tendo somente mudado os nomes dos arquivos para **nvt**.

Nos comandos anteriores (equilibração no ensaio NVT), estabilizamos a temperatura do sistema. Antes da simulação de produção, também temos de estabilizar a pressão (e, portanto, também a densidade) do sistema. A equilibração da pressão é realizada sob o ensaio NPT, em que o número de partículas, pressão e temperatura são todos constantes.
Nesta etapa de equilibração NPT utilizaremos, mais uma vez, o **grompp** e o **mdrun** assim como fizemos na etapa de equilíbrio NVT e na etapa de minimização de energia. Note que agora estamos incluindo a opção **-t** na linha de comando para incluir o arquivo de *checkpoint* gerado durante a simulação de equilíbrio NVT.

```
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
```

* **-f**: indicamos o arquivo contendo os parâmetros a serem computados e aplicados pelo programa **genion** no próximo comando;
* **-p**: indicamos o arquivo de topologia a ser utilizado e modificado;
* **-o**: arquivo de saída que será utilizado no próximo passo contendo as configurações aplicadas nesta etapa
* **-c** e **-r**: são os arquivos contendo a estrutura do sistema tal com as suas restrições. Nas versões mais novas do GROMACS é necessário adicionar esses dois argumentos;
* **-t**: arquivo de checkpoint para continuarmos os resultados gerados pela etapa de NVT.

```
gmx mdrun -v -deffnm npt
```

* **-v**: *verbose mode*, escreve na tela os passos da minimização;
* **-s**: arquivo pré compilado com informações de coordenadas, topologia e parâmetros da simulação. Esse arquivo é adquirido na etapa anterior com o **grompp**;
* **-deffnm**: coloca os nomes dos arquivos de saída com o nomes indicado.

## Dinâmica Molecular

Após a conclusão das duas fases de equilíbrio, o sistema está equilibrado, tanto em relação a temperatura quanto a pressão desejada. Agora estamos prontos para liberar as restrições de posições e executar a etapa de simulação de DM de produção para coleta dos dados a serem analisados. 
O processo é como vimos antes, vamos fazer uso do arquivo de checkpoint (que neste caso agora contém as informações de pressão, além das informações de temperatura) do “grompp” e do “mdrun”.
 
 ```
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr 
```

* **-f**: indicamos o arquivo contendo os parâmetros a serem computados e aplicados pelo programa **genion** no próximo comando;
* **-c**: arquivo com a configuração atual do sistema, originado do passo anterior;
* **-p**: indicamos o arquivo de topologia a ser utilizado e modificado;
* **-o**: arquivo de saída que será utilizado no próximo passo contendo as configurações aplicadas nesta etapa
* **-c** e **-r**: são os arquivos contendo a estrutura do sistema tal com as suas restrições. Nas versões mais novas do GROMACS é necessário adicionar esses dois argumentos;
* **-t**: arquivo de checkpoint para continuarmos os resultados gerados pela etapa de NVT.

```
gmx mdrun -v -deffnm md 
```

* **-v**: *verbose mode*, escreve na tela os passos da minimização;
* **-s**: arquivo pré compilado com informações de coordenadas, topologia e parâmetros da simulação. Esse arquivo é adquirido na etapa anterior com o **grompp**;
* **-deffnm**: coloca os nomes dos arquivos de saída com o nomes indicado.

Ao final da simulação, são gerados quatro arquivos importante: a trajetória (**.trr** e **.xtc**), o qual contém as coordenadas de todos os átomos ao longo do tempo (velocidades); a estrutura que guarda as coordenadas e velocidades do último passo (**.gro**); a energia, salva valores de energia, pressão, entre outros (**.edr**); arquivo de checkpoint caso seja necessário mais tempo de simulação (**.cpt**).
