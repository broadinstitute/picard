/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.arrays.illumina;

public class InfiniumTransformation {

    private int version;
    private float offsetX;
    private float offsetY;
    private float scaleX;
    private float scaleY;
    private float shear;
    private float theta;
    private float reserved1;
    private float reserved2;
    private float reserved3;
    private float reserved4;
    private float reserved5;
    private float reserved6;

    public int getVersion() {
        return version;
    }

    public void setVersion(int version) {
        this.version = version;
    }

    public float getOffsetX() {
        return offsetX;
    }

    public void setOffsetX(float offsetX) {
        this.offsetX = offsetX;
    }

    public float getOffsetY() {
        return offsetY;
    }

    public void setOffsetY(float offsetY) {
        this.offsetY = offsetY;
    }

    public float getScaleX() {
        return scaleX;
    }

    public void setScaleX(float scaleX) {
        this.scaleX = scaleX;
    }

    public float getScaleY() {
        return scaleY;
    }

    public void setScaleY(float scaleY) {
        this.scaleY = scaleY;
    }

    public float getShear() {
        return shear;
    }

    public void setShear(float shear) {
        this.shear = shear;
    }

    public float getTheta() {
        return theta;
    }

    public void setTheta(float theta) {
        this.theta = theta;
    }

    public float getReserved1() {
        return reserved1;
    }

    public void setReserved1(float reserved1) {
        this.reserved1 = reserved1;
    }

    public float getReserved2() {
        return reserved2;
    }

    public void setReserved2(float reserved2) {
        this.reserved2 = reserved2;
    }

    public float getReserved3() {
        return reserved3;
    }

    public void setReserved3(float reserved3) {
        this.reserved3 = reserved3;
    }

    public float getReserved4() {
        return reserved4;
    }

    public void setReserved4(float reserved4) {
        this.reserved4 = reserved4;
    }

    public float getReserved5() {
        return reserved5;
    }

    public void setReserved5(float reserved5) {
        this.reserved5 = reserved5;
    }

    public float getReserved6() {
        return reserved6;
    }

    public void setReserved6(float reserved6) {
        this.reserved6 = reserved6;
    }
}
