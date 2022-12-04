/**
 * @name Projekt 2 - Jednoducha shlukova analyza: 2D nejblizsi soused
 * @author Zdeněk Němec <xnemec0d@stud.fit.vutbr.cz>
 */

/**
 * Kostra programu pro 2. projekt IZP 2022/23
 *
 * Jednoducha shlukova analyza: 2D nejblizsi soused.
 * Single linkage
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h> // sqrtf
#include <limits.h> // INT_MAX

/*****************************************************************
 * Ladici makra. Vypnout jejich efekt lze definici makra
 * NDEBUG, napr.:
 *   a) pri prekladu argumentem prekladaci -DNDEBUG
 *   b) v souboru (na radek pred #include <assert.h>
 *      #define NDEBUG
 */
#ifdef NDEBUG
#define debug(s)
#define dfmt(s, ...)
#define dint(i)
#define dfloat(f)
#else

// vypise ladici retezec
#define debug(s) printf("- %s\n", s)

// vypise formatovany ladici vystup - pouziti podobne jako printf
#define dfmt(s, ...) printf(" - "__FILE__":%u: "s"\n",__LINE__,__VA_ARGS__)

// vypise ladici informaci o promenne - pouziti dint(identifikator_promenne)
#define dint(i) printf(" - " __FILE__ ":%u: " #i " = %d\n", __LINE__, i)

// vypise ladici informaci o promenne typu float - pouziti
// dfloat(identifikator_promenne)
#define dfloat(f) printf(" - " __FILE__ ":%u: " #f " = %g\n", __LINE__, f)

#endif

/*****************************************************************
 * Deklarace potrebnych datovych typu:
 *
 * TYTO DEKLARACE NEMENTE
 *
 *   struct obj_t - struktura objektu: identifikator a souradnice
 *   struct cluster_t - shluk objektu:
 *      pocet objektu ve shluku,
 *      kapacita shluku (pocet objektu, pro ktere je rezervovano
 *          misto v poli),
 *      ukazatel na pole shluku.
 */

typedef struct obj_t {
    int id;
    float x;
    float y;
} obj_t;

typedef struct cluster_t {
    int size;
    int capacity;
    struct obj_t *obj;
} cluster_t;

/*****************************************************************
 * Deklarace potrebnych funkci.
 *
 * PROTOTYPY FUNKCI NEMENTE
 *
 * IMPLEMENTUJTE POUZE FUNKCE NA MISTECH OZNACENYCH 'TODO'
 *
 */

/*
 Inicializace shluku 'c'. Alokuje pamet pro cap objektu (kapacitu).
 Ukazatel NULL u pole objektu znamena kapacitu 0.
*/
void init_cluster(struct cluster_t *c, int cap)
{
    assert(c != NULL);
    assert(cap >= 0);

    c->size = 0;
    c->capacity = cap;
    c->obj = malloc(sizeof(obj_t) * cap);

    if (c == NULL) c->capacity = 0;

}

/*
 Odstraneni vsech objektu shluku a inicializace na prazdny shluk.
 */
void clear_cluster(struct cluster_t *c)
{
    assert(c);
    free(c->obj);
    init_cluster(c, 0);
}

/// Chunk of cluster objects. Value recommended for reallocation.
const int CLUSTER_CHUNK = 10;

/*
 Zmena kapacity shluku 'c' na kapacitu 'new_cap'.
 */
struct cluster_t *resize_cluster(struct cluster_t *c, int new_cap)
{
    // TUTO FUNKCI NEMENTE
    assert(c);
    assert(c->capacity >= 0);
    assert(new_cap >= 0);

    if (c->capacity >= new_cap)
        return c;

    size_t size = sizeof(struct obj_t) * new_cap;

    void *arr = realloc(c->obj, size);
    if (arr == NULL)
        return NULL;

    c->obj = (struct obj_t*)arr;
    c->capacity = new_cap;
    return c;
}

/*
 Prida objekt 'obj' na konec shluku 'c'. Rozsiri shluk, pokud se do nej objekt
 nevejde.
 */
void append_cluster(struct cluster_t *c, struct obj_t obj)
{
    assert(c);

    if (c->size < c->capacity) {
        c->obj[c->size] = obj;
        c->size += 1;
    }
    else {
        c = resize_cluster(c, (c->capacity)+CLUSTER_CHUNK);
        c->obj[c->size] = obj;
        c->size += 1;
    }
}

/*
 Seradi objekty ve shluku 'c' vzestupne podle jejich identifikacniho cisla.
 */
void sort_cluster(struct cluster_t *c);

/*
 Do shluku 'c1' prida objekty 'c2'. Shluk 'c1' bude v pripade nutnosti rozsiren.
 Objekty ve shluku 'c1' budou serazeny vzestupne podle identifikacniho cisla.
 Shluk 'c2' bude nezmenen.
 */
void merge_clusters(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c2 != NULL);

    for (int i = 0; i < c2->size; i++) {
        append_cluster(c1, c2->obj[i]);
    }
    
    sort_cluster(c1);
}

/**********************************************************************/
/* Prace s polem shluku */

/*
 Odstrani shluk z pole shluku 'carr'. Pole shluku obsahuje 'narr' polozek
 (shluku). Shluk pro odstraneni se nachazi na indexu 'idx'. Funkce vraci novy
 pocet shluku v poli.
*/
int remove_cluster(struct cluster_t *carr, int narr, int idx)
{
    assert(idx < narr);
    assert(narr > 0);

    for (int i = idx; i < narr-1; i++) {
        cluster_t tmp = carr[i];
        carr[i] = carr[i+1];
        carr[i+1] = tmp;
    }

    carr = realloc(carr, (sizeof(cluster_t)*(narr-1)));
    
    return narr-1;

}

/*
 Pocita Euklidovskou vzdalenost mezi dvema objekty.
 */
float obj_distance(struct obj_t *o1, struct obj_t *o2)
{
    assert(o1 != NULL);
    assert(o2 != NULL);

    float x = (o1->x - o2->x);
    x *= x;

    float y = (o1->y - o2->y);
    y *= y;

    float euclidean_distance = sqrtf(x + y);

    return euclidean_distance;
}

/*
 Pocita vzdalenost dvou shluku.
*/
float cluster_distance(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c1->size > 0);
    assert(c2 != NULL);
    assert(c2->size > 0);

    float euclidean_distance = obj_distance(&(c1->obj[0]), &(c2->obj[0]));
    float euclidean_distance_tmp;

    for (int j = 0; j < c1->size; j++) {

        for (int k = 0; k < c2->size; k++) {
            euclidean_distance_tmp = obj_distance(&(c1->obj[j]), &(c2->obj[k]));

            if (euclidean_distance_tmp < euclidean_distance) {
                euclidean_distance = euclidean_distance_tmp;
            }
        }
    }

    return euclidean_distance;
}

/*
 Funkce najde dva nejblizsi shluky. V poli shluku 'carr' o velikosti 'narr'
 hleda dva nejblizsi shluky. Nalezene shluky identifikuje jejich indexy v poli
 'carr'. Funkce nalezene shluky (indexy do pole 'carr') uklada do pameti na
 adresu 'c1' resp. 'c2'.
*/
void find_neighbours(struct cluster_t *carr, int narr, int *c1, int *c2)
{
    assert(narr > 0);

    float distance = cluster_distance(&(carr[0]), &(carr[1]));
    *c1 = 0;
    *c2 = 1;

    for (int j = 0; j < narr; j++) {
        for (int k = j+1; k < narr; k++) {
            if (cluster_distance(&(carr[j]), &(carr[k])) < distance) {
                distance = cluster_distance(&(carr[j]), &(carr[k]));
                *c1 = j;
                *c2 = k;
            }
        }
    }
}

// pomocna funkce pro razeni shluku
static int obj_sort_compar(const void *a, const void *b)
{
    // TUTO FUNKCI NEMENTE
    const struct obj_t *o1 = (const struct obj_t *)a;
    const struct obj_t *o2 = (const struct obj_t *)b;
    if (o1->id < o2->id) return -1;
    if (o1->id > o2->id) return 1;
    return 0;
}

/*
 Razeni objektu ve shluku vzestupne podle jejich identifikatoru.
*/
void sort_cluster(struct cluster_t *c)
{
    // TUTO FUNKCI NEMENTE
    qsort(c->obj, c->size, sizeof(struct obj_t), &obj_sort_compar);
}

/*
 Tisk shluku 'c' na stdout.
*/
void print_cluster(struct cluster_t *c)
{
    // TUTO FUNKCI NEMENTE
    for (int i = 0; i < c->size; i++)
    {
        if (i) putchar(' ');
        printf("%d[%g,%g]", c->obj[i].id, c->obj[i].x, c->obj[i].y);
    }
    putchar('\n');
}

/*
 Ze souboru 'filename' nacte objekty. Pro kazdy objekt vytvori shluk a ulozi
 jej do pole shluku. Alokuje prostor pro pole vsech shluku a ukazatel na prvni
 polozku pole (ukalazatel na prvni shluk v alokovanem poli) ulozi do pameti,
 kam se odkazuje parametr 'arr'. Funkce vraci pocet nactenych objektu (shluku).
 V pripade nejake chyby uklada do pameti, kam se odkazuje 'arr', hodnotu NULL.
*/
int load_clusters(char *filename, struct cluster_t **arr)
{
    assert(arr != NULL);

    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Was not able to open file");
        return -1;
    }

    char trash1[6];
    char trash2[2];
    float countt;

    if (fscanf(f, "%5s %c %f", trash1, trash2, &countt) != 3) {
        fprintf(stderr, "first line in file is in incorrect format\n");
        return -1;

    }
    int count = countt;
    if (count != countt) {
        fprintf(stderr, "Count is not int\n");
        return -1;
    }

    if (count < 1) {
        fprintf(stderr, "there are no objects or object count is negative\n");
        return -1;
    }

    *arr = malloc(sizeof(cluster_t)*count);

    int x;
    int y;
    for (int i = 0; i < count; i++) {
        obj_t obj;
        float id;
        int idd;

        if (fscanf(f, "%f %f %f", &id, &obj.x, &obj.y) != 3) {
            fprintf(stderr, "id or coordinates are not number\n");
            return -1;
        }

        idd = id;
        if (idd != id) {
            fprintf(stderr, "ID is not int\n");
            return -1;
        }

        if ((obj.x > 1000) || (obj.y > 1000) || (obj.x < 0) || (obj.y < 0)) {
            fprintf(stderr, "Coordinations out of range\n");
            return -1;
        }

        x = obj.x;
        y = obj.y;

        if ((x != obj.x) || (y != obj.y)) {
            fprintf(stderr, "Floating point number has been found in coordinations\n");
            return -1;
        }

        obj.id = idd;

        if ((obj.id > INT_MAX) || (obj.id < INT_MIN)) {
            fprintf(stderr, "ID out of INT MIN/MAX range\n");
            return -1;
        }

        init_cluster(*arr+i, CLUSTER_CHUNK);
        append_cluster(*arr+i, obj);

        if (*arr+i == NULL) *arr = NULL;
    }

    fclose(f);

    if (*arr != NULL) return count;
    return 0;
}

/*
 Tisk pole shluku. Parametr 'carr' je ukazatel na prvni polozku (shluk).
 Tiskne se prvnich 'narr' shluku.
*/
void print_clusters(struct cluster_t *carr, int narr)
{
    printf("Clusters:\n");
    for (int i = 0; i < narr; i++)
    {
        printf("cluster %d: ", i);
        print_cluster(&(carr[i]));
    }
}

int check_ids_duplicity(cluster_t *carr, int narr) {
    for (int j = 0; j < narr; j++) {
        for (int k = j+1; k < narr; k++) {
            if (carr[j].obj[0].id == carr[k].obj[0].id) {
                return 1;
            }
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    if(!((argc>1)&&(argc<4))) {
        fprintf(stderr, "1 or 2 parameters expected, got %d\n", argc-1);
        return 1;
        }

    char *final_count_char;
    float final_count_float;
    int final_count;

    if (argc == 2) {final_count = 1;}
    else {
        final_count_char = argv[2];
        final_count_float = atof(final_count_char);
        final_count = atoi(final_count_char);
        if (final_count_float != final_count) {
            fprintf(stderr, "parameter n is not int\n");
            return 1;
        }
        }
    if (final_count < 1) {
        fprintf(stderr, "parameter is not number or is lower than 1\n");
    }

    cluster_t *clusters;

    int count = load_clusters(argv[1], &clusters);
    if (count == -1) return 1;

    if (check_ids_duplicity(clusters, count) == 1) {
        fprintf(stderr, "duplicit IDs has been found\n");
        return 1;
    }

    int *index_one = malloc(sizeof(int));
    int *index_two = malloc(sizeof(int));

    while (count != final_count)
    {

        find_neighbours(clusters, count, index_one, index_two);

        merge_clusters(&(clusters[*index_one]), &(clusters[*index_two]));

        count = remove_cluster(clusters, count, *index_two);

    }
    
    print_clusters(clusters, count);
    free(index_one);
    free(index_two);
    free(clusters);
    return 0;
}