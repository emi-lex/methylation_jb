def compare_chr(chr1, p1, chr2, p2):
    return chr1 < chr2 or (chr1 == chr2 and p1 < p2)

def binary_search(chromosomes_arr, positions_arr, c, p, lg="le"):
    if lg == "le":
        left = 0
        right = len(chromosomes_arr)
        while right - left > 1:
            mid = (right + left) // 2
            if compare_chr(c, p, chromosomes_arr[mid], positions_arr[mid]):
                right = mid
            else:
                left = mid
        return left
    if lg == "ge":
        left = -1
        right = len(chromosomes_arr) - 1
        while right - left > 1:
            mid = (right + left) // 2
            if compare_chr(chromosomes_arr[mid], positions_arr[mid], c, p):
                left = mid
            else:
                right = mid
        return right