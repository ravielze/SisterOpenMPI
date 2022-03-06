## Penjelasan Program
Garis besar kerja program adalah master node akan menerima inputan terlebih dahulu (kernel matriks dan input matriks). Lalu selanjutnya master node akan membagi input matriks ke masing-masing node sejumlah jumlah_input_matriks/jumlah_node dengan openMPI. Setelah masing-masing node mendapatkan "beban kerja"-nya mereka akan mengerjakan kerjaannya dengan skema paralelisasi threading dengan openMP. Setelah masing-masing node menyelesaikan kerjaannya, mereka akan melakukan sorting masing-masing dan melemparkan hasil sorting tersebut ke "partner"-nya. Nantinya partner-partner ini akan menggabungkan array-nya sendiri dan array yang diterima, selanjutnya array hasil gabungan akan dilemparkan ke partner selanjutnya. Hal ini terus dilakukan hingga semua array terkumpul di node 0 (master).

## Analisa Perbandingan Runtime

## Perbedaan Hasil Program Paralel dan Program Serial

## Variasi Parameter Node dan Thread
1. Node: 2, Thread: 5

2. Node: 2, Thread: 16

3. Node: 3, Thread: 5

4. Node: 3, Thread: 16

5. Node: 4, Thread: 5

6. Node: 4, Thread: 16

## Analisa Waku Eksekusi

## Author:
13519002 - Steven Nataniel Kodyat
13519007 - Muhammad Tito Prakasa
13519044 - Kinantan Arya Bagaspati