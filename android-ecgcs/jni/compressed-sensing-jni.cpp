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
 *
 */
#include <string.h>
#include <jni.h>
#include "compressed-sensing-inl.h"

/* This is a trivial JNI example where we use a native method
 * to return a new VM String. See the corresponding Java source
 * file located at:
 *
 *   apps/samples/hello-jni/project/src/com/example/HelloJni/HelloJni.java
 */

extern "C" {

jlong
Java_com_nesl_ecgcs_EcgCs_encodeBernoulli(JNIEnv* env, jobject obj, 
jlong lfsr_state, jfloatArray x_java, jfloatArray y_java) {
	int N = env->GetArrayLength(x_java);
	int K = env->GetArrayLength(y_java);

	float* x = env->GetFloatArrayElements(x_java, 0);
	float* y = env->GetFloatArrayElements(y_java, 0);

	unsigned int new_lfsr_state = encode_bernoulli(lfsr_state, N, K, x, y);

	env->ReleaseFloatArrayElements(y_java, y, 0);
	env->ReleaseFloatArrayElements(x_java, x, 0);
  return new_lfsr_state;
}

jlong
Java_com_nesl_ecgcs_EcgCs_decodeDwt(JNIEnv* env, jobject obj,
jlong lfsr_state, jfloatArray y_java, jfloatArray xhat_java) {
	int K = env->GetArrayLength(y_java);
	int N = env->GetArrayLength(xhat_java);

	float* y = env->GetFloatArrayElements(y_java, 0);
	float* xhat = env->GetFloatArrayElements(xhat_java, 0);

	int new_lfsr_state = decode_dwt(lfsr_state, K, N, y, xhat);

	env->ReleaseFloatArrayElements(y_java, y, 0);
	env->ReleaseFloatArrayElements(xhat_java, xhat, 0);
  return new_lfsr_state;
}



}  // end extern "C"
