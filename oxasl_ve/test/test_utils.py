import numpy as np

from oxasl_ve import two_to_mac, mac_to_two, veslocs_to_enc

TWO1 = np.array([
    [0, 0, 0, 0],
    [0, 1, 0, 0],
    [90, 2, -23, 22],
    [90, 3, -23, 22],
    [0, 2, -11.5, -26.5],
    [0, 3, -11.5, -26.5],
    [20.95578, 2, -17.56445, -3.19896],
    [20.95578, 3, -17.56445, -3.19896],
], dtype=np.float)

MAC1 = np.array([
    [-0.5, 44.5, 0, 0, -3.71299, 1.42479],
    [0, 0, -19.0, -4, -9.69503, 3.72029],
    [0, 0, 270, 270, 69.04422, 69.04422],
    [22.5, 22.5, 7.5, 7.5, 7.18275, 7.18275],
], dtype=np.float)

MAC2 = np.array([
    [0, 20, 0, 0, 20, 20],
    [0, 0, 0, 20, 20, 20],
    [0, 0, 270, 270, 135, 135],
    [10, 10, 10, 10, 4.472136, 4.472136],
], dtype=np.float)

def test_two_to_mac():
    mac, imlist = two_to_mac(TWO1)
    assert(np.all(np.isclose(MAC1, mac, atol=1e-4)))
    assert(np.all(imlist == np.arange(8) - 1))

def test_mac_to_two():
    two, imlist = mac_to_two(MAC1)
    assert(np.all(np.isclose(TWO1, two, atol=1e-4)))
    assert(np.all(imlist == np.arange(8) - 1))

def test_roundtrip():
    two, imlist = mac_to_two(MAC2)
    mac, imlist2 = two_to_mac(two)
    print(MAC2)
    print(two)
    print(mac)
    assert(np.all(np.isclose(mac, MAC2, atol=1e-4)))
    assert(np.all(imlist == imlist2))

def test_roundtrip2():
    two, imlist = mac_to_two(MAC1)
    mac, imlist2 = two_to_mac(two)
    assert(np.all(np.isclose(mac, MAC1, atol=1e-4)))
    assert(np.all(imlist == imlist2))
