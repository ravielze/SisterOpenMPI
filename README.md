## Penjelasan Program

Garis besar kerja program adalah master node akan menerima inputan terlebih dahulu (kernel matriks dan input matriks). Kemudian program berjalan denan urutan pengerjaan sesuai yang diberikan spesifikasi tugas, dengan no 2 dan 3 digabung:

1. Master node akan membagi input matriks ke masing-masing node sejumlah jumlah_input_matriks/jumlah_node dengan openMPI. Pengiriman dan penerimaan menggunakan MPI_Send. Hal yang dikirim ialah integer yang bernilai banyaknya matrix yang dikirim ke tiap node, serta masing masing isi matrix yang akan diproses tiap node.

2. Setelah masing-masing node mendapatkan "beban kerja"-nya mereka akan mengerjakan pekerjaannya dengan skema paralelisasi threading dengan openMP. Pada step ini kami menggabungkan konvolusi dengan mencari selisih cell terbesar dan terkecil tiap matrix, karena kami coba lebih cepat dari apabila dipisah. Idenya adalah setiap menghitung hasil konvolusi tepat 1 cell pada matrix keluaran, ubah nilai max dan min dari matrix tersebut. Hal inilah yang juga memotivasi kami untuk mengadakan skema paralelisasi di atas proses konvolusi, karena proses ini mengakses cell yang sama sehingga kami nilai tidak akan menghasilkan speedup signifikan bila proses tersebut dipecah. Kami sajikan 3 kode berbeda yang berbeda syntax namun esensi mirip, yakni dilakukan paralelisasi pada for (jmlMatrix * row * col) karena jika pada salah satu aspek saja ditakutkan terjadi ketidakseimbangan pembagian beban pada salah satu dari 5 atau 16 thread (semakin besar aspek yang dibagi akan semakin seimbang jika dibagi 5 atau 16). Tak lupa juga untuk menambahkan reduction max dan min pada array maxArray dan minArray.

3. Setelah masing-masing node menyelesaikan kerjaannya, mereka akan melakukan sorting masing-masing dan melemparkan hasil sorting tersebut ke "partner"-nya. Sorting yang dilakukan pada masing-masing node kami nilai legal karena jika ingin dipararelkan lagi maka sudah bukan sorting dengan openMPI lagi kategorinya, karena harus menggunakan openMP. Nantinya partner-partner ini akan menggabungkan array-nya sendiri dan array yang diterima, menggunakan fungsi merge_array, selanjutnya array hasil gabungan akan dilemparkan ke partner selanjutnya. Penentuan node mana yang melempar/menerima ialah dari banyaknya bit 0 di akhir representasi biner id node ini. Hal ini terus dilakukan hingga semua array terkumpul di node 0 (master).

## Analisa Perbandingan Runtime

Terdapat 5 test case (termasuk TC pada docs tubes) yang diberikan pada spesifikasi tugas. TC 0 sampai 4 memiliki ukuran kasar (estimasi) berturut-turut (2, 3, 3), (5, 33, 36), (15, 100, 100), (100, 666, 50), (24, 5000, 50), degan masing masing triplet menyatakan ukuran kernel, banyak input, dan ukuran matrix input. Cek kompleksitas sederhana, apabila 3 komponen triplet di atas dimisalkan (x,y,z), maka proram ini memiliki kompleksitas O(x*x*y*z*z + y log(y)). Pada eksekusi program di perangkat lokal dengan spesifikasi (
    ISI MAS TITO
), kami mendapati bahwa program paralel kami lebih cepat baru pada eksekusi TC2 dengan 3 node dan 5 thread. Sebetulnya hal ini masih dalam ekspektasi kami, karena pada TC yang kecil, cost yang dibutuhkan untuk membagi dan komunikasi secara umum antar node, serta inisiasi thread, tidak sepadan dengan keuntungan yang didapat dari pekerjaan yang dilakukan 2 sampai 4 kali lipat lebih cepat hanya pada sebagian kecil dari proses penghitungan. Speedup yang dirasakan semakin signifikan seiring naiknya kompleksitas dan ukuran TC. Penghitungan singkat kompleksitas program apabila dijalankan dengan n node dan t thread adalah O(xxyzz/nt + ylog(y)/n), tentunya mengabaikan konstanta "cold start" threading dan komunikasi antar thread.

## Perbedaan Hasil Program Paralel dan Program Serial

Dari segi hasil, kami mendapati program paralel kami menghasilkan nilai yang sama dengan program serial. Meskipun demikian, kami pernah mengalami perbedaan, tentunya apabila tidak diberikan beberapa batasan pada parameter pragma jika di dalamnya terdapat variabel yang dapat diakses lebih dari satu thread sekaligus. Sebelum solusi yang sekarang diberikan reduction, karena blok kode max/min terdiri dari 2 operasi yakni compare dan assign, bisa saja saat dijalankan 2 thread sekaligus akan terjadi perbedaan hasil akhir yang tertimpa karena tidak atomik. Misalkan var shared = 1, thread 1 sedang menghitung max(shared,3), sementara thread 2 max(shared, 2), ketika mulai, kedua thread akan masuk ke assigning karena 2,3 > 1, namun bila thread 2 selesai lebih lama, shared akan bernilai 2 sehingga salah. Reduksi yang dilakukan program kita hanya pada max/min karena penjumlahan konvolusi tidak kami paralelkan (berjalan pada maximal 1 thread). Proses sorting juga tidak mungkin menghasilkan perbedaan hasil karena setiap saatnya tidak ada variabel yang dapat diakses 2 node sekaligus.

## Variasi Parameter Node dan Thread
| TC | Node | Thread | Time |
|----|------|--------|------|
|  1 |   2  |    5   | kekw |
|  1 |   2  |   16   | kekw |
|  1 |   3  |    5   | kekw |
|  1 |   3  |   16   | kekw |
|  1 |   4  |    5   | kekw |
|  1 |   4  |   16   | kekw |
|  2 |   2  |    5   | kekw |
|  2 |   2  |   16   | kekw |
|  2 |   3  |    5   | kekw |
|  2 |   3  |   16   | kekw |
|  2 |   4  |    5   | kekw |
|  2 |   4  |   16   | kekw |
|  3 |   2  |    5   | kekw |
|  3 |   2  |   16   | kekw |
|  3 |   3  |    5   | kekw |
|  3 |   3  |   16   | kekw |
|  3 |   4  |    5   | kekw |
|  3 |   4  |   16   | kekw |
|  4 |   2  |    5   | kekw |
|  4 |   2  |   16   | kekw |
|  4 |   3  |    5   | kekw |
|  4 |   3  |   16   | kekw |
|  4 |   4  |    5   | kekw |
|  4 |   4  |   16   | kekw |

## Analisa Waktu Eksekusi

## Author:
13519002 - Steven Nataniel Kodyat
13519007 - Muhammad Tito Prakasa
13519044 - Kinantan Arya Bagaspati