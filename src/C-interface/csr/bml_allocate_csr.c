#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_allocate_csr.h"
#include "bml_types_csr.h"

#include <stdio.h>
#include <math.h>

/** allocate hash table.
 *
 * \ingroup allocate_group
 *
 * \param tsize - the initial hash table size.
 */
csr_row_index_hash_t *
csr_noinit_table(
    const int tsize)
{
    int i, lwr;
    const int alloc_size =
        INIT_SLOT_STORAGE_SIZE >= tsize ? INIT_SLOT_STORAGE_SIZE : tsize;
    static int powers[] =
        { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
    static int powersof2[] =
        { 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768,
        65536, 131072, 262144
    };
    for (i = 1; lwr = (powersof2[i] * 2) / 3, lwr < alloc_size; i++);
    const int space = powersof2[i];
    const int space_minus1 = space - 1;

  /** create table object */
    csr_row_index_hash_t *table =
        bml_noinit_allocate_memory(sizeof(csr_row_index_hash_t));

  /** allocate the space for member variables */
    table->size_ = 0;
    table->space_ = space;
    table->capacity_ = space;
    table->space_minus1_ = space_minus1;
    table->slot_storage_space_ = table->capacity_;
  /** allocate array of slot pointers and initialize slot pointers to NULL */
    table->Slots_ =
        bml_noinit_allocate_memory(sizeof(csr_hash_slot_t *) * space);

    for (i = 0; i < space; i++)
    {
        table->Slots_[i] = NULL;
    }

  /** Allocate memory for storing data */
    csr_hash_slot_t *newstorage =
        bml_noinit_allocate_memory(sizeof(csr_hash_slot_t) *
                                   table->slot_storage_space_);

    table->slot_storage_ = newstorage;
  /** initialize slot_ptr_ */
    table->slot_ptr_ = &newstorage[0];

    return table;
}

/** insert key into hash table.
 *
 * \ingroup allocate_group
 *
 * \param table - the hash table.
 * \param key - key to be inserted
 */
void
csr_table_insert(
    csr_row_index_hash_t * table,
    const int key)
{
    const int size = table->size_;
    //reallocate storage if needed (not used yet - DOK)
/*
    if ((size & table->space_minus1_) == 0)
    {
        table->slot_storage_space_ += table->capacity_;
        csr_hash_slot_t *newstorage =
            bml_reallocate_memory(table->slot_storage_,
                                  sizeof(csr_hash_slot_t) *
                                  table->slot_storage_space_);

        table->slot_storage_ = newstorage;
        table->slot_ptr_ = &newstorage[size];
    }
*/
    csr_hash_slot_t *slot_ptr = table->slot_ptr_;

    slot_ptr->key = key;
    slot_ptr->value = size;

    const int index = (int) hash_key_index(key, table->space_minus1_);
    slot_ptr->link = table->Slots_[index];
    table->Slots_[index] = slot_ptr;
    table->size_++;
    table->slot_ptr_++;
}

/** Get the corresponding value for a given key.
 *
 * \ingroup allocate_group
 *
 * \param table - the hash table.
 * \param key - key to be inserted
 */
void *
csr_table_lookup(
    csr_row_index_hash_t * table,
    const int key)
{
    const int index = (int) hash_key_index(key, table->space_minus1_);

    struct csr_hash_slot_t *const slot = table->Slots_[index];
    struct csr_hash_slot_t *p;
    for (p = slot; p; p = p->link)
    {
        if (p->key == key)
        {
            return &p->value;
        }
    }
    return NULL;
}

/** Deallocate hash table.
 *
 * \ingroup allocate_group
 *
 * \param table - the hash table.
 */
void
csr_deallocate_table(
    csr_row_index_hash_t * table)
{
    /** delete allocated slots */
    bml_free_memory(table->slot_storage_);
    bml_free_memory(table->Slots_);
    bml_free_memory(table);
}

/** Deallocate csr row.
 *
 * \ingroup allocate_group
 *
 * \param row - the csr row.
 */
void
csr_deallocate_row(
    csr_sparse_row_t * row)
{
    bml_free_memory(row->vals_);
    bml_free_memory(row->cols_);
    bml_free_memory(row);
}

/** Deallocate a matrix.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void
bml_deallocate_csr(
    bml_matrix_csr_t * A)
{
//    csr_deallocate_table(A->table_);
    /** deallocate row data */
    const int n = A->N_;
    for (int i = 0; i < n; i++)
    {
        csr_deallocate_row((A->data_)[i]);
    }
    bml_free_memory(A->data_);
//    bml_free_memory(A->lvarsgid_);
//    bml_deallocate_domain(A->domain);
//    bml_deallocate_domain(A->domain2);
    bml_free_memory(A);
}

/** Clear a matrix.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void
bml_clear_csr(
    bml_matrix_csr_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_clear_csr_single_real(A);
            break;
        case double_real:
            bml_clear_csr_double_real(A);
            break;
        case single_complex:
            bml_clear_csr_single_complex(A);
            break;
        case double_complex:
            bml_clear_csr_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Allocate the csr matrix.
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param M The maximum number of non-zeros per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_csr_t *
bml_noinit_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_csr_t *A = NULL;

    switch (matrix_precision)
    {
        case single_real:
            A = bml_noinit_matrix_csr_single_real(matrix_dimension,
                                                  distrib_mode);
            break;
        case double_real:
            A = bml_noinit_matrix_csr_double_real(matrix_dimension,
                                                  distrib_mode);
            break;
        case single_complex:
            A = bml_noinit_matrix_csr_single_complex(matrix_dimension,
                                                     distrib_mode);
            break;
        case double_complex:
            A = bml_noinit_matrix_csr_double_complex(matrix_dimension,
                                                     distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return A;
}

/** Allocate the zero matrix. (Currently assumes sequential case only)
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_csr_t *
bml_zero_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_csr_t *A = NULL;

    switch (matrix_precision)
    {
        case single_real:
            A = bml_zero_matrix_csr_single_real(N, M, distrib_mode);
            break;
        case double_real:
            A = bml_zero_matrix_csr_double_real(N, M, distrib_mode);
            break;
        case single_complex:
            A = bml_zero_matrix_csr_single_complex(N, M, distrib_mode);
            break;
        case double_complex:
            A = bml_zero_matrix_csr_double_complex(N, M, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return A;
}

/** Allocate a banded random matrix. (Currently assumes sequential case only)
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */

bml_matrix_csr_t *
bml_banded_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_banded_matrix_csr_single_real(N, M, distrib_mode);
            break;
        case double_real:
            return bml_banded_matrix_csr_double_real(N, M, distrib_mode);
            break;
        case single_complex:
            return bml_banded_matrix_csr_single_complex(N, M, distrib_mode);
            break;
        case double_complex:
            return bml_banded_matrix_csr_double_complex(N, M, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}

/** Allocate a random matrix. (Currently assumes sequential case only)
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */

bml_matrix_csr_t *
bml_random_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_random_matrix_csr_single_real(N, M, distrib_mode);
            break;
        case double_real:
            return bml_random_matrix_csr_double_real(N, M, distrib_mode);
            break;
        case single_complex:
            return bml_random_matrix_csr_single_complex(N, M, distrib_mode);
            break;
        case double_complex:
            return bml_random_matrix_csr_double_complex(N, M, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}

/** Allocate the identity matrix. (Currently assumes sequential case only)
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_csr_t *
bml_identity_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_identity_matrix_csr_single_real(N, M, distrib_mode);
            break;
        case double_real:
            return bml_identity_matrix_csr_double_real(N, M, distrib_mode);
            break;
        case single_complex:
            return bml_identity_matrix_csr_single_complex(N, M, distrib_mode);
            break;
        case double_complex:
            return bml_identity_matrix_csr_double_complex(N, M, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}

/** Update the csr matrix domain.
 *
 * \ingroup allocate_group
 *
 * \param A Matrix with domain
 * \param localPartMin first part on each rank
 * \param localPartMax last part on each rank
 * \param nnodesInPart number of nodes per part
 */
void
bml_update_domain_csr(
    bml_matrix_csr_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart)
{
    LOG_ERROR("bml_update_domain_csr not implemented\n");
}
