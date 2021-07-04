#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED_10ai {

static constexpr std::array<sunindextype, 5103> dxdotdw_rowvals_SPARCED_10ai_ = {
    0, 0, 1, 3, 4, 37, 38, 39, 42, 42, 42, 43, 43, 43, 47, 47, 50, 54, 53, 55, 56, 60, 62, 66, 69, 72, 75, 78, 83, 84, 84, 87, 89, 89, 92, 94, 94, 97, 101, 103, 106, 109, 109, 109, 111, 124, 127, 131, 134, 141, 142, 143, 147, 155, 156, 162, 164, 166, 167, 169, 170, 157, 171, 171, 158, 643, 644, 645, 648, 648, 649, 649, 649, 649, 650, 650, 650, 650, 653, 659, 659, 662, 662, 662, 663, 666, 668, 670, 670, 673, 676, 676, 676, 676, 679, 681, 683, 685, 687, 692, 694, 694, 697, 700, 701, 704, 706, 707, 710, 710, 710, 710, 713, 715, 715, 718, 720, 722, 723, 723, 725, 727, 34, 35, 172, 172, 159, 159, 761, 642, 642, 160, 160, 173, 28, 27, 29, 174, 161, 1, 3, 4, 27, 28, 29, 34, 35, 37, 38, 39, 42, 43, 47, 50, 53, 54, 55, 56, 60, 62, 66, 69, 72, 75, 78, 83, 84, 87, 89, 92, 94, 97, 101, 103, 106, 109, 111, 124, 127, 131, 134, 141, 142, 143, 147, 155, 156, 157, 158, 159, 160, 161, 162, 164, 166, 167, 169, 170, 171, 172, 173, 174, 642, 643, 644, 645, 648, 649, 650, 653, 659, 662, 663, 666, 668, 670, 673, 676, 679, 681, 683, 685, 687, 692, 694, 697, 700, 701, 704, 706, 707, 710, 713, 715, 718, 720, 722, 723, 725, 727, 761, 761, 771, 772, 761, 771, 772, 771, 771, 772, 1, 1, 1, 2, 1, 2, 1, 2, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 3, 33, 2, 35, 36, 2, 35, 36, 30, 31, 31, 30, 31, 39, 39, 42, 76, 39, 42, 76, 39, 40, 39, 40, 39, 40, 41, 40, 41, 40, 42, 77, 40, 42, 77, 40, 41, 42, 42, 63, 42, 63, 42, 63, 43, 43, 44, 43, 44, 43, 43, 44, 45, 44, 45, 45, 46, 60, 45, 46, 60, 47, 47, 48, 47, 48, 47, 48, 49, 48, 49, 49, 51, 60, 49, 51, 60, 50, 50, 50, 54, 52, 54, 52, 54, 54, 52, 56, 56, 57, 56, 57, 56, 56, 57, 58, 57, 58, 58, 59, 60, 58, 59, 60, 60, 64, 60, 64, 64, 64, 61, 62, 61, 62, 61, 62, 62, 53, 53, 65, 53, 65, 53, 65, 66, 66, 67, 66, 67, 66, 66, 67, 68, 67, 68, 60, 68, 74, 60, 68, 74, 69, 69, 70, 69, 70, 69, 70, 55, 55, 71, 55, 71, 55, 71, 72, 72, 73, 72, 73, 72, 73, 75, 75, 47, 60, 60, 60, 45, 78, 79, 45, 78, 79, 49, 78, 80, 49, 78, 80, 58, 78, 81, 58, 78, 81, 68, 78, 82, 68, 78, 82, 78, 83, 84, 85, 83, 84, 85, 85, 86, 86, 87, 88, 86, 87, 88, 86, 89, 90, 86, 89, 90, 86, 90, 91, 91, 92, 93, 91, 92, 93, 91, 94, 95, 91, 94, 95, 91, 95, 96, 96, 97, 98, 96, 97, 98, 96, 98, 99, 89, 99, 100, 89, 99, 100, 91, 99, 100, 96, 101, 102, 96, 101, 102, 101, 102, 140, 96, 103, 104, 96, 103, 104, 96, 104, 105, 91, 106, 107, 91, 106, 107, 91, 107, 108, 108, 109, 110, 108, 109, 110, 108, 111, 112, 108, 111, 112, 108, 112, 113, 113, 114, 113, 114, 114, 115, 116, 114, 115, 116, 114, 117, 114, 117, 115, 117, 118, 115, 117, 118, 117, 119, 117, 119, 115, 119, 120, 115, 119, 120, 119, 121, 122, 119, 121, 122, 119, 122, 123, 123, 124, 125, 123, 124, 125, 123, 125, 126, 123, 127, 128, 123, 127, 128, 123, 128, 129, 126, 130, 126, 130, 130, 131, 132, 130, 131, 132, 130, 132, 133, 133, 134, 135, 133, 134, 135, 94, 135, 136, 94, 135, 136, 96, 135, 136, 129, 137, 129, 137, 101, 135, 138, 101, 135, 138, 101, 137, 139, 101, 137, 139, 109, 141, 144, 109, 141, 144, 109, 142, 145, 109, 142, 145, 109, 143, 146, 109, 143, 146, 111, 147, 148, 111, 147, 148, 113, 147, 148, 109, 147, 149, 109, 147, 149, 121, 123, 127, 129, 124, 126, 109, 115, 109, 115, 89, 91, 155, 162, 188, 155, 162, 188, 164, 188, 194, 164, 188, 194, 162, 188, 196, 162, 188, 196, 155, 196, 214, 155, 196, 214, 188, 214, 188, 214, 166, 188, 197, 166, 188, 197, 156, 197, 215, 156, 197, 215, 188, 189, 215, 188, 189, 215, 169, 188, 195, 169, 188, 195, 167, 188, 198, 167, 188, 198, 156, 198, 216, 156, 198, 216, 188, 190, 216, 188, 190, 216, 156, 166, 189, 156, 166, 189, 164, 189, 199, 164, 189, 199, 162, 189, 200, 162, 189, 200, 155, 200, 215, 155, 200, 215, 166, 189, 201, 166, 189, 201, 156, 201, 217, 156, 201, 217, 189, 217, 189, 217, 169, 189, 202, 169, 189, 202, 167, 189, 203, 167, 189, 203, 156, 203, 218, 156, 203, 218, 189, 190, 218, 189, 190, 218, 156, 167, 190, 156, 167, 190, 164, 190, 204, 164, 190, 204, 162, 190, 206, 162, 190, 206, 155, 206, 216, 155, 206, 216, 166, 190, 207, 166, 190, 207, 156, 207, 218, 156, 207, 218, 169, 190, 205, 169, 190, 205, 167, 190, 208, 167, 190, 208, 156, 208, 219, 156, 208, 219, 190, 219, 190, 219, 162, 225, 717, 162, 163, 164, 226, 717, 164, 165, 167, 227, 717, 167, 168, 157, 170, 191, 157, 170, 191, 170, 191, 209, 170, 191, 209, 157, 209, 220, 157, 209, 220, 191, 220, 191, 220, 158, 171, 192, 158, 171, 192, 192, 221, 192, 221, 158, 210, 221, 171, 192, 210, 159, 172, 193, 159, 172, 193, 193, 222, 193, 222, 159, 211, 222, 159, 211, 222, 172, 193, 211, 172, 193, 211, 173, 186, 173, 186, 174, 187, 174, 187, 162, 225, 717, 163, 225, 717, 164, 226, 717, 165, 226, 717, 167, 227, 717, 168, 227, 717, 162, 175, 162, 175, 162, 164, 176, 162, 164, 176, 162, 166, 177, 162, 166, 177, 162, 167, 178, 162, 167, 178, 164, 179, 164, 179, 164, 166, 180, 164, 166, 180, 164, 167, 181, 164, 167, 181, 166, 167, 182, 166, 167, 182, 167, 183, 167, 183, 155, 175, 196, 155, 175, 196, 155, 176, 194, 155, 176, 194, 155, 177, 197, 155, 177, 197, 156, 177, 200, 156, 177, 200, 155, 178, 198, 155, 178, 198, 156, 178, 206, 156, 178, 206, 156, 180, 199, 156, 180, 199, 156, 181, 204, 156, 181, 204, 156, 182, 207, 156, 182, 207, 156, 183, 208, 156, 183, 208, 170, 184, 170, 184, 157, 184, 209, 157, 184, 209, 172, 185, 172, 185, 159, 185, 211, 159, 185, 211, 160, 186, 212, 160, 186, 212, 160, 212, 223, 160, 212, 223, 161, 187, 213, 161, 187, 213, 161, 213, 224, 161, 213, 224, 194, 228, 195, 229, 196, 230, 214, 231, 197, 232, 215, 233, 198, 234, 216, 235, 199, 236, 202, 237, 200, 238, 203, 239, 218, 240, 204, 241, 205, 242, 206, 243, 207, 244, 208, 245, 219, 246, 209, 247, 220, 248, 221, 249, 210, 250, 222, 251, 211, 252, 212, 253, 213, 254, 223, 255, 224, 256, 194, 228, 195, 229, 196, 230, 214, 231, 197, 232, 215, 233, 198, 234, 216, 235, 199, 236, 202, 237, 200, 238, 203, 239, 218, 240, 204, 241, 205, 242, 206, 243, 207, 244, 208, 245, 219, 246, 209, 247, 220, 248, 221, 249, 210, 250, 222, 251, 211, 252, 212, 253, 213, 254, 223, 255, 224, 256, 194, 261, 643, 195, 262, 643, 196, 263, 643, 214, 264, 643, 197, 265, 643, 215, 266, 643, 198, 267, 643, 216, 268, 643, 199, 269, 643, 202, 270, 643, 200, 271, 643, 203, 272, 643, 218, 273, 643, 204, 274, 643, 205, 275, 643, 206, 276, 643, 207, 277, 643, 208, 278, 643, 219, 279, 643, 209, 280, 643, 220, 281, 643, 221, 282, 643, 210, 283, 643, 222, 284, 643, 211, 285, 643, 212, 286, 643, 213, 287, 643, 223, 288, 643, 224, 289, 643, 194, 261, 643, 195, 262, 643, 196, 263, 643, 214, 264, 643, 197, 265, 643, 215, 266, 643, 198, 267, 643, 216, 268, 643, 199, 269, 643, 202, 270, 643, 200, 271, 643, 203, 272, 643, 218, 273, 643, 204, 274, 643, 205, 275, 643, 206, 276, 643, 207, 277, 643, 208, 278, 643, 219, 279, 643, 209, 280, 643, 220, 281, 643, 221, 282, 643, 210, 283, 643, 222, 284, 643, 211, 285, 643, 212, 286, 643, 213, 287, 643, 223, 288, 643, 224, 289, 643, 228, 319, 229, 320, 230, 321, 231, 322, 232, 323, 233, 324, 234, 325, 235, 326, 236, 327, 237, 328, 238, 329, 239, 330, 240, 331, 241, 332, 242, 333, 243, 334, 244, 335, 245, 336, 246, 337, 247, 338, 248, 339, 249, 340, 250, 341, 251, 342, 252, 343, 253, 344, 254, 345, 255, 346, 256, 347, 290, 319, 291, 320, 292, 321, 293, 322, 294, 323, 295, 324, 296, 325, 297, 326, 298, 327, 299, 328, 300, 329, 301, 330, 302, 331, 303, 332, 304, 333, 305, 334, 306, 335, 307, 336, 308, 337, 309, 338, 310, 339, 311, 340, 312, 341, 313, 342, 314, 343, 315, 344, 316, 345, 317, 346, 318, 347, 290, 319, 291, 320, 292, 321, 293, 322, 294, 323, 295, 324, 296, 325, 297, 326, 298, 327, 299, 328, 300, 329, 301, 330, 302, 331, 303, 332, 304, 333, 305, 334, 306, 335, 307, 336, 308, 337, 309, 338, 310, 339, 311, 340, 312, 341, 313, 342, 314, 343, 315, 344, 316, 345, 317, 346, 318, 347, 194, 290, 195, 291, 196, 292, 214, 293, 197, 294, 215, 295, 198, 296, 216, 297, 199, 298, 202, 299, 200, 300, 203, 301, 218, 302, 204, 303, 205, 304, 206, 305, 207, 306, 208, 307, 219, 308, 209, 309, 220, 310, 221, 311, 210, 312, 222, 313, 211, 314, 212, 315, 213, 316, 223, 317, 224, 318, 319, 381, 646, 320, 382, 646, 321, 383, 646, 322, 384, 646, 323, 385, 646, 324, 386, 646, 325, 387, 646, 326, 388, 646, 327, 389, 646, 328, 390, 646, 329, 391, 646, 330, 392, 646, 331, 393, 646, 332, 394, 646, 333, 395, 646, 334, 396, 646, 335, 397, 646, 336, 398, 646, 337, 399, 646, 338, 400, 646, 339, 401, 646, 340, 402, 646, 341, 403, 646, 342, 404, 646, 343, 405, 646, 348, 406, 646, 349, 407, 646, 350, 408, 646, 351, 409, 646, 319, 381, 646, 320, 382, 646, 321, 383, 646, 322, 384, 646, 323, 385, 646, 324, 386, 646, 325, 387, 646, 326, 388, 646, 327, 389, 646, 328, 390, 646, 329, 391, 646, 330, 392, 646, 331, 393, 646, 332, 394, 646, 333, 395, 646, 334, 396, 646, 335, 397, 646, 336, 398, 646, 337, 399, 646, 338, 400, 646, 339, 401, 646, 340, 402, 646, 341, 403, 646, 342, 404, 646, 343, 405, 646, 348, 406, 646, 349, 407, 646, 350, 408, 646, 351, 409, 646, 228, 352, 646, 229, 353, 646, 230, 354, 646, 231, 355, 646, 232, 356, 646, 233, 357, 646, 234, 358, 646, 235, 359, 646, 236, 360, 646, 237, 361, 646, 238, 362, 646, 239, 363, 646, 240, 364, 646, 241, 365, 646, 242, 366, 646, 243, 367, 646, 244, 368, 646, 245, 369, 646, 246, 370, 646, 247, 371, 646, 248, 372, 646, 249, 373, 646, 250, 374, 646, 251, 375, 646, 252, 376, 646, 257, 377, 646, 258, 378, 646, 259, 379, 646, 260, 380, 646, 228, 352, 646, 229, 353, 646, 230, 354, 646, 231, 355, 646, 232, 356, 646, 233, 357, 646, 234, 358, 646, 235, 359, 646, 236, 360, 646, 237, 361, 646, 238, 362, 646, 239, 363, 646, 240, 364, 646, 241, 365, 646, 242, 366, 646, 243, 367, 646, 244, 368, 646, 245, 369, 646, 246, 370, 646, 247, 371, 646, 248, 372, 646, 249, 373, 646, 250, 374, 646, 251, 375, 646, 252, 376, 646, 257, 377, 646, 258, 378, 646, 259, 379, 646, 260, 380, 646, 228, 410, 648, 229, 411, 648, 230, 412, 648, 231, 413, 648, 232, 414, 648, 233, 415, 648, 234, 416, 648, 235, 417, 648, 236, 418, 648, 237, 419, 648, 238, 420, 648, 239, 421, 648, 240, 422, 648, 241, 423, 648, 242, 424, 648, 243, 425, 648, 244, 426, 648, 245, 427, 648, 246, 428, 648, 247, 429, 648, 248, 430, 648, 249, 431, 648, 250, 432, 648, 251, 433, 648, 252, 434, 648, 257, 435, 648, 258, 436, 648, 259, 437, 648, 260, 438, 648, 228, 410, 648, 229, 411, 648, 230, 412, 648, 231, 413, 648, 232, 414, 648, 233, 415, 648, 234, 416, 648, 235, 417, 648, 236, 418, 648, 237, 419, 648, 238, 420, 648, 239, 421, 648, 240, 422, 648, 241, 423, 648, 242, 424, 648, 243, 425, 648, 244, 426, 648, 245, 427, 648, 246, 428, 648, 247, 429, 648, 248, 430, 648, 249, 431, 648, 250, 432, 648, 251, 433, 648, 252, 434, 648, 257, 435, 648, 258, 436, 648, 259, 437, 648, 260, 438, 648, 228, 439, 651, 229, 440, 651, 230, 441, 651, 231, 442, 651, 232, 443, 651, 233, 444, 651, 234, 445, 651, 235, 446, 651, 236, 447, 651, 237, 448, 651, 238, 449, 651, 239, 450, 651, 240, 451, 651, 241, 452, 651, 242, 453, 651, 243, 454, 651, 244, 455, 651, 245, 456, 651, 246, 457, 651, 247, 458, 651, 248, 459, 651, 249, 460, 651, 250, 461, 651, 251, 462, 651, 252, 463, 651, 257, 464, 651, 258, 465, 651, 259, 466, 651, 260, 467, 651, 228, 439, 651, 229, 440, 651, 230, 441, 651, 231, 442, 651, 232, 443, 651, 233, 444, 651, 234, 445, 651, 235, 446, 651, 236, 447, 651, 237, 448, 651, 238, 449, 651, 239, 450, 651, 240, 451, 651, 241, 452, 651, 242, 453, 651, 243, 454, 651, 244, 455, 651, 245, 456, 651, 246, 457, 651, 247, 458, 651, 248, 459, 651, 249, 460, 651, 250, 461, 651, 251, 462, 651, 252, 463, 651, 257, 464, 651, 258, 465, 651, 259, 466, 651, 260, 467, 651, 228, 468, 653, 229, 469, 653, 230, 470, 653, 231, 471, 653, 232, 472, 653, 233, 473, 653, 234, 474, 653, 235, 475, 653, 236, 476, 653, 237, 477, 653, 238, 478, 653, 239, 479, 653, 240, 480, 653, 241, 481, 653, 242, 482, 653, 243, 483, 653, 244, 484, 653, 245, 485, 653, 246, 486, 653, 247, 487, 653, 248, 488, 653, 249, 489, 653, 250, 490, 653, 251, 491, 653, 252, 492, 653, 257, 493, 653, 258, 494, 653, 259, 495, 653, 260, 496, 653, 228, 468, 653, 229, 469, 653, 230, 470, 653, 231, 471, 653, 232, 472, 653, 233, 473, 653, 234, 474, 653, 235, 475, 653, 236, 476, 653, 237, 477, 653, 238, 478, 653, 239, 479, 653, 240, 480, 653, 241, 481, 653, 242, 482, 653, 243, 483, 653, 244, 484, 653, 245, 485, 653, 246, 486, 653, 247, 487, 653, 248, 488, 653, 249, 489, 653, 250, 490, 653, 251, 491, 653, 252, 492, 653, 257, 493, 653, 258, 494, 653, 259, 495, 653, 260, 496, 653, 381, 497, 662, 382, 498, 662, 383, 499, 662, 384, 500, 662, 385, 501, 662, 386, 502, 662, 387, 503, 662, 388, 504, 662, 389, 505, 662, 390, 506, 662, 391, 507, 662, 392, 508, 662, 393, 509, 662, 394, 510, 662, 395, 511, 662, 396, 512, 662, 397, 513, 662, 398, 514, 662, 399, 515, 662, 400, 516, 662, 401, 517, 662, 402, 518, 662, 403, 519, 662, 404, 520, 662, 405, 521, 662, 406, 522, 662, 407, 523, 662, 408, 524, 662, 409, 525, 662, 381, 497, 662, 382, 498, 662, 383, 499, 662, 384, 500, 662, 385, 501, 662, 386, 502, 662, 387, 503, 662, 388, 504, 662, 389, 505, 662, 390, 506, 662, 391, 507, 662, 392, 508, 662, 393, 509, 662, 394, 510, 662, 395, 511, 662, 396, 512, 662, 397, 513, 662, 398, 514, 662, 399, 515, 662, 400, 516, 662, 401, 517, 662, 402, 518, 662, 403, 519, 662, 404, 520, 662, 405, 521, 662, 406, 522, 662, 407, 523, 662, 408, 524, 662, 409, 525, 662, 381, 497, 661, 382, 498, 661, 383, 499, 661, 384, 500, 661, 385, 501, 661, 386, 502, 661, 387, 503, 661, 388, 504, 661, 389, 505, 661, 390, 506, 661, 391, 507, 661, 392, 508, 661, 393, 509, 661, 394, 510, 661, 395, 511, 661, 396, 512, 661, 397, 513, 661, 398, 514, 661, 399, 515, 661, 400, 516, 661, 401, 517, 661, 402, 518, 661, 403, 519, 661, 404, 520, 661, 405, 521, 661, 406, 522, 661, 407, 523, 661, 408, 524, 661, 409, 525, 661, 352, 526, 662, 353, 527, 662, 354, 528, 662, 355, 529, 662, 356, 530, 662, 357, 531, 662, 358, 532, 662, 359, 533, 662, 360, 534, 662, 361, 535, 662, 362, 536, 662, 363, 537, 662, 364, 538, 662, 365, 539, 662, 366, 540, 662, 367, 541, 662, 368, 542, 662, 369, 543, 662, 370, 544, 662, 371, 545, 662, 372, 546, 662, 373, 547, 662, 374, 548, 662, 375, 549, 662, 376, 550, 662, 377, 551, 662, 378, 552, 662, 379, 553, 662, 380, 554, 662, 352, 526, 662, 353, 527, 662, 354, 528, 662, 355, 529, 662, 356, 530, 662, 357, 531, 662, 358, 532, 662, 359, 533, 662, 360, 534, 662, 361, 535, 662, 362, 536, 662, 363, 537, 662, 364, 538, 662, 365, 539, 662, 366, 540, 662, 367, 541, 662, 368, 542, 662, 369, 543, 662, 370, 544, 662, 371, 545, 662, 372, 546, 662, 373, 547, 662, 374, 548, 662, 375, 549, 662, 376, 550, 662, 377, 551, 662, 378, 552, 662, 379, 553, 662, 380, 554, 662, 352, 526, 661, 353, 527, 661, 354, 528, 661, 355, 529, 661, 356, 530, 661, 357, 531, 661, 358, 532, 661, 359, 533, 661, 360, 534, 661, 361, 535, 661, 362, 536, 661, 363, 537, 661, 364, 538, 661, 365, 539, 661, 366, 540, 661, 367, 541, 661, 368, 542, 661, 369, 543, 661, 370, 544, 661, 371, 545, 661, 372, 546, 661, 373, 547, 661, 374, 548, 661, 375, 549, 661, 376, 550, 661, 377, 551, 661, 378, 552, 661, 379, 553, 661, 380, 554, 661, 410, 555, 690, 411, 556, 690, 412, 557, 690, 413, 558, 690, 414, 559, 690, 415, 560, 690, 416, 561, 690, 417, 562, 690, 418, 563, 690, 419, 564, 690, 420, 565, 690, 421, 566, 690, 422, 567, 690, 423, 568, 690, 424, 569, 690, 425, 570, 690, 426, 571, 690, 427, 572, 690, 428, 573, 690, 429, 574, 690, 430, 575, 690, 431, 576, 690, 432, 577, 690, 433, 578, 690, 434, 579, 690, 435, 580, 690, 436, 581, 690, 437, 582, 690, 438, 583, 690, 410, 555, 690, 411, 556, 690, 412, 557, 690, 413, 558, 690, 414, 559, 690, 415, 560, 690, 416, 561, 690, 417, 562, 690, 418, 563, 690, 419, 564, 690, 420, 565, 690, 421, 566, 690, 422, 567, 690, 423, 568, 690, 424, 569, 690, 425, 570, 690, 426, 571, 690, 427, 572, 690, 428, 573, 690, 429, 574, 690, 430, 575, 690, 431, 576, 690, 432, 577, 690, 433, 578, 690, 434, 579, 690, 435, 580, 690, 436, 581, 690, 437, 582, 690, 438, 583, 690, 410, 555, 658, 689, 411, 556, 658, 689, 412, 557, 658, 689, 413, 558, 658, 689, 414, 559, 658, 689, 415, 560, 658, 689, 416, 561, 658, 689, 417, 562, 658, 689, 418, 563, 658, 689, 419, 564, 658, 689, 420, 565, 658, 689, 421, 566, 658, 689, 422, 567, 658, 689, 423, 568, 658, 689, 424, 569, 658, 689, 425, 570, 658, 689, 426, 571, 658, 689, 427, 572, 658, 689, 428, 573, 658, 689, 429, 574, 658, 689, 430, 575, 658, 689, 431, 576, 658, 689, 432, 577, 658, 689, 433, 578, 658, 689, 434, 579, 658, 689, 435, 580, 658, 689, 436, 581, 658, 689, 437, 582, 658, 689, 438, 583, 658, 689, 439, 584, 690, 440, 585, 690, 441, 586, 690, 442, 587, 690, 443, 588, 690, 444, 589, 690, 445, 590, 690, 446, 591, 690, 447, 592, 690, 448, 593, 690, 449, 594, 690, 450, 595, 690, 451, 596, 690, 452, 597, 690, 453, 598, 690, 454, 599, 690, 455, 600, 690, 456, 601, 690, 457, 602, 690, 458, 603, 690, 459, 604, 690, 460, 605, 690, 461, 606, 690, 462, 607, 690, 463, 608, 690, 464, 609, 690, 465, 610, 690, 466, 611, 690, 467, 612, 690, 439, 584, 690, 440, 585, 690, 441, 586, 690, 442, 587, 690, 443, 588, 690, 444, 589, 690, 445, 590, 690, 446, 591, 690, 447, 592, 690, 448, 593, 690, 449, 594, 690, 450, 595, 690, 451, 596, 690, 452, 597, 690, 453, 598, 690, 454, 599, 690, 455, 600, 690, 456, 601, 690, 457, 602, 690, 458, 603, 690, 459, 604, 690, 460, 605, 690, 461, 606, 690, 462, 607, 690, 463, 608, 690, 464, 609, 690, 465, 610, 690, 466, 611, 690, 467, 612, 690, 439, 584, 691, 440, 585, 691, 441, 586, 691, 442, 587, 691, 443, 588, 691, 444, 589, 691, 445, 590, 691, 446, 591, 691, 447, 592, 691, 448, 593, 691, 449, 594, 691, 450, 595, 691, 451, 596, 691, 452, 597, 691, 453, 598, 691, 454, 599, 691, 455, 600, 691, 456, 601, 691, 457, 602, 691, 458, 603, 691, 459, 604, 691, 460, 605, 691, 461, 606, 691, 462, 607, 691, 463, 608, 691, 464, 609, 691, 465, 610, 691, 466, 611, 691, 467, 612, 691, 468, 613, 656, 469, 614, 656, 470, 615, 656, 471, 616, 656, 472, 617, 656, 473, 618, 656, 474, 619, 656, 475, 620, 656, 476, 621, 656, 477, 622, 656, 478, 623, 656, 479, 624, 656, 480, 625, 656, 481, 626, 656, 482, 627, 656, 483, 628, 656, 484, 629, 656, 485, 630, 656, 486, 631, 656, 487, 632, 656, 488, 633, 656, 489, 634, 656, 490, 635, 656, 491, 636, 656, 492, 637, 656, 493, 638, 656, 494, 639, 656, 495, 640, 656, 496, 641, 656, 468, 613, 656, 469, 614, 656, 470, 615, 656, 471, 616, 656, 472, 617, 656, 473, 618, 656, 474, 619, 656, 475, 620, 656, 476, 621, 656, 477, 622, 656, 478, 623, 656, 479, 624, 656, 480, 625, 656, 481, 626, 656, 482, 627, 656, 483, 628, 656, 484, 629, 656, 485, 630, 656, 486, 631, 656, 487, 632, 656, 488, 633, 656, 489, 634, 656, 490, 635, 656, 491, 636, 656, 492, 637, 656, 493, 638, 656, 494, 639, 656, 495, 640, 656, 496, 641, 656, 468, 613, 657, 469, 614, 657, 470, 615, 657, 471, 616, 657, 472, 617, 657, 473, 618, 657, 474, 619, 657, 475, 620, 657, 476, 621, 657, 477, 622, 657, 478, 623, 657, 479, 624, 657, 480, 625, 657, 481, 626, 657, 482, 627, 657, 483, 628, 657, 484, 629, 657, 485, 630, 657, 486, 631, 657, 487, 632, 657, 488, 633, 657, 489, 634, 657, 490, 635, 657, 491, 636, 657, 492, 637, 657, 493, 638, 657, 494, 639, 657, 495, 640, 657, 496, 641, 657, 253, 257, 642, 253, 257, 642, 254, 258, 642, 254, 258, 642, 255, 259, 642, 255, 259, 642, 256, 260, 642, 256, 260, 642, 344, 348, 642, 344, 348, 642, 345, 349, 642, 345, 349, 642, 346, 350, 642, 346, 350, 642, 347, 351, 642, 347, 351, 642, 661, 666, 667, 661, 666, 667, 661, 668, 757, 661, 668, 757, 667, 763, 667, 763, 757, 758, 757, 758, 667, 669, 757, 667, 669, 757, 670, 763, 764, 670, 763, 764, 671, 763, 764, 671, 763, 765, 671, 763, 765, 672, 763, 765, 670, 758, 759, 670, 758, 759, 671, 758, 759, 671, 758, 760, 671, 758, 760, 672, 758, 760, 669, 670, 733, 669, 670, 733, 669, 671, 733, 669, 671, 734, 669, 671, 734, 669, 672, 734, 670, 671, 671, 672, 672, 715, 735, 672, 715, 735, 672, 716, 735, 672, 716, 736, 672, 716, 736, 672, 717, 736, 675, 679, 739, 675, 679, 739, 674, 679, 739, 679, 680, 681, 682, 707, 708, 677, 678, 646, 717, 728, 646, 717, 728, 647, 717, 728, 646, 647, 666, 717, 729, 666, 717, 729, 665, 717, 729, 665, 666, 660, 662, 730, 660, 662, 730, 660, 661, 730, 661, 663, 731, 661, 663, 731, 662, 663, 731, 663, 717, 732, 663, 717, 732, 664, 717, 732, 663, 664, 673, 716, 754, 673, 716, 754, 673, 715, 754, 716, 717, 673, 717, 755, 673, 716, 755, 676, 717, 737, 676, 717, 737, 677, 717, 737, 676, 677, 675, 717, 674, 675, 675, 680, 756, 675, 680, 756, 674, 680, 756, 674, 715, 677, 678, 678, 679, 738, 678, 679, 738, 678, 680, 738, 675, 680, 739, 678, 681, 740, 678, 681, 740, 678, 682, 740, 675, 681, 741, 675, 681, 741, 675, 682, 741, 682, 683, 684, 682, 683, 684, 658, 659, 660, 658, 659, 660, 658, 710, 711, 658, 710, 711, 711, 713, 742, 711, 713, 742, 711, 712, 742, 712, 713, 666, 713, 714, 666, 713, 714, 691, 692, 743, 691, 692, 743, 690, 692, 743, 691, 693, 694, 691, 693, 694, 691, 697, 698, 691, 697, 698, 700, 701, 702, 700, 701, 702, 693, 698, 744, 693, 698, 744, 698, 699, 744, 693, 699, 699, 702, 745, 699, 702, 745, 702, 703, 745, 699, 703, 691, 696, 703, 691, 696, 703, 696, 704, 746, 696, 704, 746, 696, 705, 746, 704, 705, 687, 688, 687, 688, 686, 687, 686, 687, 696, 707, 747, 696, 707, 747, 696, 708, 747, 707, 717, 748, 707, 717, 748, 708, 717, 748, 706, 707, 709, 706, 707, 709, 654, 655, 654, 655, 654, 701, 722, 654, 701, 722, 655, 725, 750, 655, 725, 750, 655, 726, 750, 725, 726, 655, 723, 751, 655, 723, 751, 655, 724, 751, 723, 724, 696, 718, 752, 696, 718, 752, 696, 719, 752, 718, 719, 651, 655, 753, 651, 655, 753, 652, 655, 753, 651, 652, 649, 650, 651, 649, 650, 651, 656, 657, 645, 646, 727, 645, 646, 727, 725, 761, 762, 725, 761, 762, 690, 688, 715, 716, 716, 717, 674, 675, 695, 696, 694, 695, 720, 721, 718, 766, 718, 766, 690, 691, 690, 691, 661, 662, 661, 662, 672, 767, 768, 672, 767, 768, 694, 769, 770, 694, 769, 770, 3, 32, 696, 3, 32, 696, 32, 33, 696, 3, 33, 147, 150, 717, 147, 150, 717, 150, 151, 717, 147, 151, 141, 152, 696, 141, 152, 696, 152, 153, 696, 141, 154, 717, 141, 154, 717, 153, 154, 717, 141, 153, 85, 86, 88, 90, 91, 93, 95, 96, 98, 99, 100, 102, 104, 105, 107, 108, 110, 112, 113, 114, 115, 116, 117, 118, 119, 120, 126, 129, 130, 132, 133, 135, 136, 137, 138, 139, 140, 144, 145, 146, 148, 149, 150, 151, 152, 153, 154, 163, 165, 168, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 646, 647, 651, 652, 654, 655, 656, 657, 658, 660, 661, 664, 665, 667, 669, 671, 672, 674, 675, 677, 678, 680, 682, 684, 686, 688, 689, 690, 691, 693, 695, 696, 698, 699, 702, 703, 705, 708, 709, 711, 712, 714, 716, 717, 719, 721, 724, 726, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 762, 763, 764, 765, 766, 768, 770
};

void dxdotdw_rowvals_SPARCED_10ai(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_SPARCED_10ai_));
}
} // namespace model_SPARCED_10ai
} // namespace amici
