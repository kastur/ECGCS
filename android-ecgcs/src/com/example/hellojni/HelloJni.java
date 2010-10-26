/*
 * Copyright (C) 2009 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.example.hellojni;

import android.app.Activity;
import android.widget.TextView;
import android.os.Bundle;


public class HelloJni extends Activity
{
    /** Called when the activity is first created. */
    @Override
    public void onCreate(Bundle savedInstanceState)
    {
        super.onCreate(savedInstanceState);

        /* Create a TextView and set its content.
         * the text is retrieved by calling a native
         * function.
         */

        TextView  tv = new TextView(this);
        tv.setText("Compressive Sensing ECG");
        setContentView(tv);

        long lfsr_state = 854875398; 
        int N = 12;
        int K = 10;

        float[] x = new float[N];
        x[0] = 118;
        x[1] = 117;
        x[2] = 117;
        x[3] = 117;
        x[4] = 117;
        x[5] = 117;
        x[6] = 117;
        x[7] = 117;
        x[8] = 118;
        x[9] = 120;
        x[10] = 122;
        x[11] = 123;

        float[] xhat = new float[N];
        for (int i = 0; i < x.length; i++) {
        	xhat[i] = 0;
        }

        float[] y = new float[K];
        for (int i = 0; i < y.length; i++) {
        	y[i] = 0;
        }

        encode(lfsr_state, x, y);
        System.out.println("y = ");
        for (int i = 0; i < y.length; i++) {
        	System.out.println(y[i]);
        }
        System.out.println("");

        decode(lfsr_state, y, xhat);
        System.out.println("x, xhat = ");
        for (int i = 0; i < x.length; i++) {
        	System.out.println(x[i] + " " + xhat[i]);
        }
    }

    public native long encode(long i, float[] x, float[] y);
    public native long decode(long i, float[] y, float[] xhat);
    
    static {
        System.loadLibrary("hello-jni");
    }
}
