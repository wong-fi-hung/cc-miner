// Copyright (c) 2009-2010 Satoshi Nakamoto
// Copyright (c) 2009-2014 The Bitcoin Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "primitives/block.h"

#include "hash.h"
#include "tinyformat.h"
#include "utilstrencodings.h"
#include "crypto/common.h"

int32_t ASSETCHAINS_MAGIC = -497513811;

// default hash algorithm for block
uint256 (CBlockHeader::*CBlockHeader::hashFunction)() const = &CBlockHeader::GetSHA256DHash;

uint256 CBlockHeader::GetSHA256DHash() const
{
    return SerializeHash(*this);
}

uint256 CBlockHeader::GetVerusHash() const
{
    if (hashPrevBlock.IsNull())
        // always use SHA256D for genesis block
        return SerializeHash(*this);
    else
        return SerializeVerusHash(*this);
}

uint256 CBlockHeader::GetVerusV2Hash() const
{
    if (hashPrevBlock.IsNull())
    {
        // always use SHA256D for genesis block
        return SerializeHash(*this);
    }
    else
    {
        if (nVersion > 4)
        {
            return SerializeVerusHashV2b(*this);
        }
        else
        {
            return SerializeVerusHash(*this);
        }
    }
}

void CBlockHeader::SetSHA256DHash()
{
    CBlockHeader::hashFunction = &CBlockHeader::GetSHA256DHash;
}

void CBlockHeader::SetVerusHash()
{
    CBlockHeader::hashFunction = &CBlockHeader::GetVerusHash;
}

// returns false if unable to fast calculate the VerusPOSHash from the header. it can still be calculated from the block
// in that case. the only difference between this and the POS hash for the contest is that it is not divided by the value out
// this is used as a source of entropy
bool CBlockHeader::GetRawVerusPOSHash(uint256 &value, int32_t nHeight) const
{
    // if below the required height or no storage space in the solution, we can't get
    // a cached txid value to calculate the POSHash from the header
    if (!(CPOSNonce::NewNonceActive(nHeight) && IsVerusPOSBlock()))
        return false;
    
    // if we can calculate, this assumes the protocol that the POSHash calculation is:
    //    hashWriter << ASSETCHAINS_MAGIC;
    //    hashWriter << nNonce; (nNonce is:
    //                           (high 128 bits == low 128 bits of verus hash of low 128 bits of nonce)
    //                           (low 32 bits == compact PoS difficult)
    //                           (mid 96 bits == low 96 bits of HASH(pastHash, txid, voutnum)
    //                              pastHash is hash of height - 100, either PoW hash of block or PoS hash, if new PoS
    //                          )
    //    hashWriter << height;
    //    return hashWriter.GetHash();
    CVerusHashWriter hashWriter = CVerusHashWriter(SER_GETHASH, PROTOCOL_VERSION);

    hashWriter << ASSETCHAINS_MAGIC;
    hashWriter << nNonce;
    hashWriter << nHeight;
    value = hashWriter.GetHash();
    return true;
}

// depending on the height of the block and its type, this returns the POS hash or the POW hash
uint256 CBlockHeader::GetVerusEntropyHash(int32_t height) const
{
    uint256 retVal;
    // if we qualify as PoW, use PoW hash, regardless of PoS state
    if (GetRawVerusPOSHash(retVal, height))
    {
        // POS hash
        return retVal;
    }
    return GetHash();
}

uint256 CBlock::BuildMerkleTree(bool* fMutated) const
{
    /* WARNING! If you're reading this because you're learning about crypto
       and/or designing a new system that will use merkle trees, keep in mind
       that the following merkle tree algorithm has a serious flaw related to
       duplicate txids, resulting in a vulnerability (CVE-2012-2459).

       The reason is that if the number of hashes in the list at a given time
       is odd, the last one is duplicated before computing the next level (which
       is unusual in Merkle trees). This results in certain sequences of
       transactions leading to the same merkle root. For example, these two
       trees:

                    A               A
                  /  \            /   \
                B     C         B       C
               / \    |        / \     / \
              D   E   F       D   E   F   F
             / \ / \ / \     / \ / \ / \ / \
             1 2 3 4 5 6     1 2 3 4 5 6 5 6

       for transaction lists [1,2,3,4,5,6] and [1,2,3,4,5,6,5,6] (where 5 and
       6 are repeated) result in the same root hash A (because the hash of both
       of (F) and (F,F) is C).

       The vulnerability results from being able to send a block with such a
       transaction list, with the same merkle root, and the same block hash as
       the original without duplication, resulting in failed validation. If the
       receiving node proceeds to mark that block as permanently invalid
       however, it will fail to accept further unmodified (and thus potentially
       valid) versions of the same block. We defend against this by detecting
       the case where we would hash two identical hashes at the end of the list
       together, and treating that identically to the block having an invalid
       merkle root. Assuming no double-SHA256 collisions, this will detect all
       known ways of changing the transactions without affecting the merkle
       root.
    */
    vMerkleTree.clear();
    vMerkleTree.reserve(vtx.size() * 2 + 16); // Safe upper bound for the number of total nodes.
    for (std::vector<CTransaction>::const_iterator it(vtx.begin()); it != vtx.end(); ++it)
        vMerkleTree.push_back(it->GetHash());
    int j = 0;
    bool mutated = false;
    for (int nSize = vtx.size(); nSize > 1; nSize = (nSize + 1) / 2)
    {
        for (int i = 0; i < nSize; i += 2)
        {
            int i2 = std::min(i+1, nSize-1);
            if (i2 == i + 1 && i2 + 1 == nSize && vMerkleTree[j+i] == vMerkleTree[j+i2]) {
                // Two identical hashes at the end of the list at a particular level.
                mutated = true;
            }
            vMerkleTree.push_back(Hash(BEGIN(vMerkleTree[j+i]),  END(vMerkleTree[j+i]),
                                       BEGIN(vMerkleTree[j+i2]), END(vMerkleTree[j+i2])));
        }
        j += nSize;
    }
    if (fMutated) {
        *fMutated = mutated;
    }
    return (vMerkleTree.empty() ? uint256() : vMerkleTree.back());
}

std::vector<uint256> CBlock::GetMerkleBranch(int nIndex) const
{
    if (vMerkleTree.empty())
        BuildMerkleTree();
    std::vector<uint256> vMerkleBranch;
    int j = 0;
    for (int nSize = vtx.size(); nSize > 1; nSize = (nSize + 1) / 2)
    {
        int i = std::min(nIndex^1, nSize-1);
        vMerkleBranch.push_back(vMerkleTree[j+i]);
        nIndex >>= 1;
        j += nSize;
    }
    return vMerkleBranch;
}

uint256 CBlock::CheckMerkleBranch(uint256 hash, const std::vector<uint256>& vMerkleBranch, int nIndex)
{
    if (nIndex == -1)
        return uint256();
    for (std::vector<uint256>::const_iterator it(vMerkleBranch.begin()); it != vMerkleBranch.end(); ++it)
    {
        if (nIndex & 1)
            hash = Hash(BEGIN(*it), END(*it), BEGIN(hash), END(hash));
        else
            hash = Hash(BEGIN(hash), END(hash), BEGIN(*it), END(*it));
        nIndex >>= 1;
    }
    return hash;
}

std::string CBlock::ToString() const
{
    std::stringstream s;
    s << strprintf("CBlock(hash=%s, ver=%d, hashPrevBlock=%s, hashMerkleRoot=%s, hashReserved=%s, nTime=%u, nBits=%08x, nNonce=%s, vtx=%u)\n",
        GetHash().ToString(),
        nVersion,
        hashPrevBlock.ToString(),
        hashMerkleRoot.ToString(),
        hashReserved.ToString(),
        nTime, nBits, nNonce.ToString(),
        vtx.size());
    /*for (unsigned int i = 0; i < vtx.size(); i++)
    {
        s << "  " << vtx[i].ToString() << "\n";
    }*/
    s << "  vMerkleTree: ";
    for (unsigned int i = 0; i < vMerkleTree.size(); i++)
        s << " " << vMerkleTree[i].ToString();
    s << "\n";
    return s.str();
}
