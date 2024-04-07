package org.husonlab.fmhdist.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import jloda.util.FileUtils;

public class FileProducer {
    private class Consumer implements FileConsumer {
        private String key;
        private String currentLine;
        private String nextLine;
        private BufferedReader reader;
        private boolean isReady;
        private boolean hasError;
        private IOException e;

        public Consumer(String filename) throws IOException {
            this.reader = new BufferedReader(new InputStreamReader(FileUtils.getInputStreamPossiblyZIPorGZIP(filename)));
            this.currentLine = this.reader.readLine();
            this.nextLine = this.reader.readLine();
            this.key = filename;
        }

        public String getKey() {
            return this.key;
        }

        public void close() throws IOException {
            this.reader.close();
        }

        @Override
        public boolean isReady() {
            return this.isReady;
        }

        public void feedLine() {
            this.currentLine = this.nextLine;
            try {
                this.nextLine = reader.readLine();
            } catch (IOException e) {
                this.hasError = true;
                this.e = e;
            }
            this.isReady = true;
        }

        @Override
        public String getLine() throws IOException {
            if (this.hasError) {
                throw this.e;
            }
            this.isReady = false;
            return this.currentLine;
        }
    }

    // Caution: need to check thread safety
	private ConcurrentHashMap<String, Consumer> consumers;
	private ConcurrentHashMap<String, Boolean> waitingForRegister;
	private ConcurrentHashMap<String, IOException> exceptions;
    private boolean isClosed;

	public FileProducer() {
		this.consumers = new ConcurrentHashMap<>();
        this.waitingForRegister = new ConcurrentHashMap<>();
        this.exceptions = new ConcurrentHashMap<>();
	}

	public void registerFile(String filename) {
		this.waitingForRegister.put(filename, true);
	}

	public FileConsumer getFileConsumer(String filename) throws IOException {
        if (this.exceptions.containsKey(filename)) {
            throw this.exceptions.get(filename);
        }
		return consumers.get(filename);
	}

	public void run() {
        while (true) {
            for(Entry<String, Boolean> e : this.waitingForRegister.entrySet()) {
                try {
                    Consumer c = new Consumer(e.getKey());
                    this.consumers.put(e.getKey(), c);
                    this.waitingForRegister.remove(e.getKey());
                } catch (IOException ex) {
                    this.exceptions.put(e.getKey(), ex);
                }
            }
            for(Entry<String, Consumer> e : this.consumers.entrySet()) {
                if(!e.getValue().isReady()) {
                    e.getValue().feedLine();
                }
            }
            if (this.isClosed) {
                return;
            }
        }
	}

    public void closeConsumer(FileConsumer consumer) throws IOException {
        // TODO: this cast could break.
        Consumer c = (Consumer) consumer;
        String fileName = c.getKey();
        if (this.consumers.containsKey(fileName)) {
            this.consumers.remove(fileName);
            c.close();
        }
    }

    public void close() {
        this.isClosed = true;
    }
}
