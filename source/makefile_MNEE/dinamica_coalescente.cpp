#include<math.h>
#include<iostream>
#include<string.h>
#include<stdio.h>
#include<stdlib.h> 
#include<iomanip> 				
#include<fstream>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_roots.h>
#include "mt.h"
#include "species_richness.h"

// salvar as árvores para visualização?
// se verdadeiro, gera um arquivo dot pra cada simulação
#define SAVE_TREES 1

#define sign(x) (x > 0) - (x < 0)

using namespace std;

struct rede{
    int N;
    int width;
    int height;
    int *barrier;
    double pb;
    double P;
    int *pilha_act;
    int size_pilha;
    int *label;
    int *species;
    double range_dispersal;
};

// estrutura e função auxiliares para o cálculo de especiação em função da
// riqueza
struct richness_function_params {
    double S;
    tree *arvore;
};
double richness_function(double U, void *params);

void coalescence(struct rede *network, double U, int (*dispersal)(double), tree *arvore, char dot_fname[]);
void definenet(rede *network, char *landscape);
void save_tree_dotfile(char *filename, tree *arvore, int *orig, int width);
double estimate_spec_rate(int S, tree *arvore, bool verbose);

/*
 * dispersal kernels functions
 */
int uniform(double range){
    return ((int) (genrand_real1() * (range - 0.5)));
}
/* Computes random sample from a normal distribution with mean zero and std.
   dev. of `sigma`, using the Box-Muller transform method. */
int normal(double sigma){
    // this is wasteful, since we use two samples to get a single sample out
    return ((int) (sigma * sqrt(-2 * log(genrand_real3())) * cos(2*M_PI*genrand_real3())));
}
/* Computes random sample from a Laplace distribution with mean zero and std.
   dev. of `sigma` */
int laplace(double sigma){
    double u;
    u = - genrand_real2() + 0.5;
    return ((int) (-sigma/sqrt(2.) * sign(u)*log(1-2*abs(u))));
 }

int main(int ac, char **av)
{
    rede network;
    int config, conf, i, j, dispersal_type, S;
    double U;
    unsigned long seed;
    char *landscape_fname, *output_fname, dot_fname[100] = "";
    ofstream output_stream;
    int (*dispersal_kernel)(double);

    tree arvore;

    if (ac != 11)
    {
        cout  <<  "start the program like this:\n" << av[0]
            << " <height> <width> <U> <S> <config> <seed> <range> <dispersal kernel> <landscape file> <output file>\n"
            << endl;
        exit (-1);
    }

    j=0;
    //Inputs sao fornecidos na linha de comando e lidos pelo atoi. Os elementos na linha por ordem são:
    network.height = atoi (av[++j]); // "altura" da rede (# de linhas)
    network.width = atoi (av[++j]); // "largura" da rede (# de colunas)
    // TODO: ser mais esperto e pegar o tamanho da rede pelo arquivo de
    // entrada? Talvez desnecessário devido à interface com R
    U = atof (av[++j]); // taxa de especiação (mínima)
    S = atoi (av[++j]); // riqueza observada (ou 0 caso U não deva ser ajustado)
    config = atoi (av[++j]); // numero de simulacoes
    seed = atol(av[++j]); // semente do RNG (numero qualquer, convertido pra unsigned long, não-negativo)
    network.range_dispersal = atof (av[++j]); // Kernel de disperao
    dispersal_type = atoi(av[++j]); // dispersal kernel, either 0 (uniform) or 1 (normal) or 2 (laplacian)
    // It would be *very* nice if this function could be provided externally
    // from R, but that would require incorporation of the function through
    // more sophisticated means than simply the command line...
    landscape_fname = av[++j]; // data file containing landscape matrix
    output_fname = av[++j]; // output file name

    cout  << "#invocation: ";
    for (int i=0; i<ac; i++){
        cout << av[i] << " ";
    }
    cout << endl;

    network.N = network.width * network.height;
    network.pilha_act = new int[network.N+1];
    network.species = new int[network.N+1];
    network.barrier = new int[network.N+1];
    network.label = new int[network.N+1];

    if (dispersal_type == 0)
        dispersal_kernel = &uniform;
    else if (dispersal_type == 1){
        dispersal_kernel = &normal;
        /* Keeping parameters comparable: value of sigma which makes the normal
           distribution have the same variance as the uniform. Even then,
           results seem to indicate there is a big difference. */
        //network.range_dispersal /= sqrt(12.);
    }
    else if (dispersal_type == 2){
        dispersal_kernel = &laplace;
        /* Keeping parameters comparable: value of sigma which makes the normal
           distribution have the same variance as the uniform. Even then,
           results seem to indicate there is a big difference. */
        //network.range_dispersal /= sqrt(12.);
    }
    else {
        cout << "Dispersal kernel type not recognized. Must be one of "
            "0 (uniform) or 1 (normal) or 2 (Laplacian)  got '" << dispersal_type << "'." << endl;
        exit(EXIT_FAILURE);
    }

    // initialize RNG
    init_genrand((unsigned long) seed);


    if (S > 0){
        // calculamos a taxa de especiação que melhor ajusta a riqueza observada
        // TODO : isto é feito usando apenas uma rede, e não tomando uma
        //        média. Investigar melhor como o código do Rosindell fazia isso.
        arvore.data.SetSize(0);

        definenet(&network, landscape_fname);
        coalescence(&network, U, dispersal_kernel, &arvore, dot_fname);

        U = estimate_spec_rate(S, &arvore, 0);
        if (U < 0){
            cout << "Erro! Taxa de especiação mínima grande demais. A diferença "
                "de riqueza foi de " << U << endl;
            return 1;
        }
        // abre o arquivo de saída só em caso de sucesso
        output_stream.open(output_fname);
        output_stream << "# parameters: N = " << network.N << " U = " << U <<
            " S = " << S << " Range = " << network.range_dispersal <<
            " Dispersal kernel = " << dispersal_type << " seed = " << seed <<
            " Landscape datafile = " << landscape_fname << endl;
        output_stream << "# taxa de especiação estimada: " << U << endl;
        cout << "taxa de especiação estimada: " << U << endl;
    } else {
        output_stream.open(output_fname);
        output_stream << "# parameters: N = " << network.N << " U = " << U <<
            " Range = " << network.range_dispersal <<
            " Dispersal kernel = " << dispersal_type << " seed = " << seed <<
            " Landscape datafile = " << landscape_fname << endl;
    }

    conf = 0;
    while((++conf) <= config){
        if (SAVE_TREES)
            sprintf(dot_fname, "%s.%d.dot", output_fname, conf);
        // limpando data por via das dúvidas
        arvore.data.SetSize(0);

        definenet(&network, landscape_fname);
        coalescence(&network, U, dispersal_kernel, &arvore, dot_fname);

        // output full grid with identity of each point
        for (i = 1; i <= network.N; i++)
            if (network.barrier[i] == 2)
                output_stream << network.species[i] << " ";
        output_stream << endl;
    }
    output_stream.close();
    return(0);
}

/* Paisagem lida a partir de um arquivo de dados, contendo 
 * - 0's (não-habitat) ou
 * - 1's (habitat, fora da área observada) ou
 * - 2's (habitat dentro da área do plot, que queremos simular)
 * O arquivo é lido por linhas, e seu comprimento total deve ser igual a width
 * * height fornecidos como parâmetros. */
void definenet(rede *network, char *landscape)
{
    int i;
    ifstream landscape_file;

    landscape_file.open(landscape);

    for( i=0; i<network->N; i++ )
    {
        if(! (landscape_file >> network->barrier[i])){
            cout << "Erro: arquivo de dados contém " << i << " entradas, esperado " << network->N << endl;
            exit(EXIT_FAILURE);
        }
    }
    landscape_file.close();

    network->size_pilha = 0;
    for( i=0; i<network->N; i++ )
    {
        if (network->barrier[i] == 2)
            network->size_pilha++;
    }
}

// Função para coalescência
void coalescence(rede *network, double U, int (*dispersal_kernel)(double), tree *arvore, char dot_fname[])
{
    int i, indice, k, x, y, xc, yc, xl, yl, size, *pilha_act, *position,
        *orig, *pilha, position_aux;
    
    pilha = new int[network->size_pilha+1];
    pilha_act = new int[network->N+1];
    position = new int[network->N+1];
    orig = new int[network->N+1];

    /***  AQUI ****/
    /*** Define o vetor inicial com posições de 0 a J-1, que chamamos pilha. 
      A pilha ativa começa no presente com J indivíduos e termina com um indivíduo. 
      A pilha pode ser vita como uma matriz em que o primeiro elemento da linha 1 tem valor 1, o segundo 2 ...
      ou seja preenchendo uma matriz by row***/
    size = 0;
    for( i=0; i<network->N; i++ )
        if (network->barrier[i] == 2)
        {
            size++;
            pilha[size] = i;
        }

    for( i=1; i<=size; i++ )
    {
        pilha_act[i] = i;
        position[i] = pilha[i];
        orig[i] = pilha[i];
    }

    // estrutura de árvore
    arvore->minspecsetup = U;
    arvore->data.SetSize(2 * (network->size_pilha) + 1);
    arvore->enddata = 0; // 0 is reserved as null
    for (int i = 0 ; i <= network->size_pilha; i++)
    {
        // looping over the entire survey area
        // setup the data variable
        arvore->enddata ++;
        arvore->data[arvore->enddata].setup(true); // generation number = 1
    }

    /***Sorteia o elemento da pilha que vai sofrer a ação ***/
    while( size>1 )
    {
        indice = (int)(genrand_real1()*size + 1);
        if(U > genrand_real1()) /*** se houve especiação o elemento é removido***/
        {
            for( i=indice; i<size; i++ )
                pilha_act[i] = pilha_act[i+1];
            size--;
        }
        else /***Caso contrário houve dispersão ***/
        {
            arvore->data[pilha_act[indice]].inc_steps();
            /* the lines below used to have by-column notation. It is
               consistent with the data file and with C-style array to use
               by-row notation. This should also be consistent with the
               definition of position_aux below. */
            xc = position[pilha_act[indice]] % network->width; //converte o valor da pilha
            yc = position[pilha_act[indice]] / network->width;

            do {
                do {
                    x = (*dispersal_kernel)(network->range_dispersal);
                    xl = xc + x;
                } while ((xl < 0) || (xl >= network->width));
                do {
                    y = (*dispersal_kernel)(network->range_dispersal);
                    yl = yc + y;
                } while ((yl < 0) || (yl >= network->height));
                /*** Aplica condições de contorno (no caso algo 
                  com reflexiva: se o sorteio joga o individuo para fora um
                  novo sorteio é feito até ele ficar dentro***/
                // Evita que o indivíduo seja "pai" de si mesmo
                if (xc == 0 && yc == 0)
                    continue;

                // see definition of xc, yc above
                position_aux = xl + yl * network->width;
                /*** Aplica a máscara de fragmentação se o sorteio joga o indivíduo na matriz
                  um novo sorteio é feito até que caia  em área de habitat***/
            } while (network->barrier[position_aux] == 0);

            /***Registra a nova posição do indivíduo. A posição no presente está em origin ***/
            position[pilha_act[indice]] = position_aux;
            // Se a posição já está ocupada coalesce
            for( k=1; k<=size; k++ )
                if(k != indice && position[pilha_act[k]] == position[pilha_act[indice]] )
                {
                    arvore->enddata++;
                    arvore->data[arvore->enddata].setup(false);
                    arvore->data[pilha_act[indice]].set_parent(arvore->enddata);
                    arvore->data[pilha_act[k]].set_parent(arvore->enddata);
                    pilha_act[k] = arvore->enddata;
                    position[arvore->enddata] = position[pilha_act[indice]];

                    for(int j = indice; j < size; j++ )
                        pilha_act[j] = pilha_act[j+1];
                    size--;

                    break;
                }
        }
    }

    // vamos salvar as árvores?
    if (dot_fname != "")
        save_tree_dotfile(dot_fname, arvore, orig, network->width);

    for(i = 1; i <= arvore->enddata; i++)
        network->species[orig[i]] = arvore->get_greatest_parent(i);

    delete[] pilha;
    delete[] pilha_act;
    delete[] position;
    delete[] orig;
}

/*
 * Salva uma árvore (ou floresta) em um arquivo DOT.
 *
 * As "folhas" (as pontas do grafo) são contornadas em vermelho e mostram a
 * posição inicial daquele indivíduo (x, y). Os labels dos demais nós não são
 * informativos (representam a numeração interna). Os labels das arestas
 * indicam o número de passos de cada indivíduo até o evento de coalescência.
 *
 * Para obter imagens a partir disso, use o Graphviz e rode "dot -O -Tpng
 * arquivo.dot".
 */
void save_tree_dotfile(char *filename, tree *arvore, int *orig, int width)
{
    int xc, yc;
    ofstream output;

    output.open(filename);
    output << "digraph arvore {" << endl;
    for (int i = 1; i < arvore->enddata; i++){
        if (arvore->data[i].get_parent() != 0)
            output << i << "->" << arvore->data[i].get_parent() << " [label=" << arvore->data[i].get_steps() << "];" << endl;
        if (arvore->data[i].get_root()){
            xc = orig[i] % width;
            yc = orig[i] / width;
            output << i << " [label=\"" << xc << ", " << yc << "\" color=red];" << endl;
        }
    }
    output << "}" << endl;
}

// função auxiliar
double richness_function(double U, void *params){
    struct richness_function_params *p = (struct richness_function_params *) params;
    return p->arvore->get_richness(U) - p->S;
}

double estimate_spec_rate(int S, tree *arvore, bool verbose){
    // adaptado da documentação oficial da GNU Scientific Library
    // https://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;
    double x_lo = 5. * arvore->minspecsetup, x_hi = 0.999;
    gsl_function F;
    struct richness_function_params params = { (double) S, arvore };

    F.function = &richness_function;
    F.params = &params;

    // taxa de especiação não foi pequena o suficiente
    if (richness_function(x_lo, &params) > 0)
        return -richness_function(x_lo, &params);

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);

    if (verbose){
        printf("using %s method\n", gsl_root_fsolver_name (s));
        printf("%5s [%9s, %9s] %9s %10s %9s\n",
                "iter", "lower", "upper", "root", "err", "err(est)");
    }

    do {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.0001);

        if (verbose){
            if (status == GSL_SUCCESS)
                printf ("Converged:\n");
            printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
                    iter, x_lo, x_hi, r, x_hi - x_lo);
        }
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);
    return r;
}
